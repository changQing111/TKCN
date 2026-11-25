#!/usr/bin/env python3
"""Fast read-only mmap querying of panel 31-mers against sample KMC databases.

Replaces building a panel KMC database and running `kmc_tools simple intersect`
+ `kmc_tools transform dump` by invoking `query_mmap` directly on sample reads
KMC databases and emitting TKCN-compatible headerless `*.realcount.tsv` tables.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


COMPLEMENT = str.maketrans("ACGT", "TGCA")
FALLBACK_QUERY_MMAP = (
    Path.home()
    / "my_data/Project/TaFTL_4A_junction31_population385_20261022/query_mmap"
)


def canonical_kmer(seq):
    seq = seq.strip().upper()
    if len(seq) != 31 or set(seq) - set("ACGT"):
        raise ValueError("Expected an ACGT 31-mer, got: {!r}".format(seq))
    rc = seq.translate(COMPLEMENT)[::-1]
    return seq if seq <= rc else rc


def normalize_db_prefix(raw_prefix):
    text = str(raw_prefix)
    for ext in (".kmc_pre", ".kmc_suf"):
        if text.endswith(ext):
            text = text[: -len(ext)]
            break
    return Path(text)


def infer_sample_name(db_prefix):
    name = normalize_db_prefix(db_prefix).name
    for suffix in (".reads.k31", ".k31"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if not name:
        raise ValueError("Cannot infer sample name from DB prefix: {}".format(db_prefix))
    return name


def validate_safe_component(value, label):
    if not value or "/" in value or "\\" in value or value in (".", ".."):
        raise ValueError("Invalid {} for path component: {!r}".format(label, value))


def resolve_query_mmap(explicit_path=None):
    candidates = []
    if explicit_path:
        candidates.append(Path(explicit_path))
    else:
        env_bin = os.environ.get("QUERY_MMAP_BIN")
        if env_bin:
            candidates.append(Path(env_bin))
        repo_bin = Path(__file__).resolve().parent / "query_mmap" / "query_mmap"
        candidates.append(repo_bin)
        which_bin = shutil.which("query_mmap")
        if which_bin:
            candidates.append(Path(which_bin))
        candidates.append(FALLBACK_QUERY_MMAP)

    for candidate in candidates:
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return candidate

    checked = ", ".join(str(p) for p in candidates)
    raise FileNotFoundError(
        "query_mmap executable not found (checked: {}). "
        "Build it with `bash src/query_mmap/build.sh` or set QUERY_MMAP_BIN / --query-mmap.".format(
            checked
        )
    )


def check_kmc_db(db_prefix):
    pre = Path(str(db_prefix) + ".kmc_pre")
    suf = Path(str(db_prefix) + ".kmc_suf")
    if not pre.is_file() or pre.stat().st_size == 0:
        raise FileNotFoundError("Missing or empty KMC prefix file: {}".format(pre))
    if not suf.is_file() or suf.stat().st_size == 0:
        raise FileNotFoundError("Missing or empty KMC suffix file: {}".format(suf))


def is_fasta_file(path):
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return stripped.startswith(">")
    raise ValueError("Panel file is empty: {}".format(path))


def prepare_panel_fasta(panel_path, tmp_dir):
    """Return a FASTA path suitable for query_mmap.

    FASTA inputs are used directly without copying. TSV inputs (panel TSV with
    `kmer`/`sequence` header or headerless 1/2-column TSV) are converted to a
    deduplicated canonical FASTA in `tmp_dir`.
    """
    panel_path = Path(panel_path)
    if not panel_path.is_file():
        raise FileNotFoundError("Panel file not found: {}".format(panel_path))
    if is_fasta_file(panel_path):
        return panel_path

    kmers = set()
    with panel_path.open(newline="") as handle:
        lines = [
            line.strip()
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not lines:
        raise ValueError("Panel TSV is empty: {}".format(panel_path))

    first_cols = lines[0].split("\t")
    col_index = 0
    start_idx = 0
    lowered = [c.strip().lower() for c in first_cols]
    if "kmer" in lowered:
        col_index = lowered.index("kmer")
        start_idx = 1
    elif "sequence" in lowered:
        col_index = lowered.index("sequence")
        start_idx = 1

    for line_no, line in enumerate(lines[start_idx:], start=start_idx + 1):
        cols = line.split("\t")
        if col_index >= len(cols):
            raise ValueError(
                "Panel TSV {} line {} missing k-mer column {}".format(
                    panel_path, line_no, col_index + 1
                )
            )
        kmers.add(canonical_kmer(cols[col_index]))

    if not kmers:
        raise ValueError("No k-mers found in panel TSV: {}".format(panel_path))

    fasta_path = Path(tmp_dir) / "panel_queries.k31.fa"
    with fasta_path.open("w") as out:
        for idx, seq in enumerate(sorted(kmers), start=1):
            out.write(">panel_{:07d}\n{}\n".format(idx, seq))
    return fasta_path


def load_membership(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError("Membership TSV not found: {}".format(path))
    gene_kmers = defaultdict(set)
    gene_roles = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"kmer", "gene_id"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                "Membership TSV {} must contain columns: {}".format(
                    path, ", ".join(sorted(required))
                )
            )
        has_role = "role" in (reader.fieldnames or [])
        row_count = 0
        for row in reader:
            row_count += 1
            gene_id = row["gene_id"].strip()
            validate_safe_component(gene_id, "gene_id")
            kmer = canonical_kmer(row["kmer"])
            role = row["role"].strip() if has_role and row.get("role") else "panel"
            validate_safe_component(role, "role")
            if gene_id in gene_roles and gene_roles[gene_id] != role:
                raise ValueError(
                    "Conflicting roles for gene_id {} in {}".format(gene_id, path)
                )
            gene_roles[gene_id] = role
            gene_kmers[gene_id].add(kmer)
    if row_count == 0:
        raise ValueError("Membership TSV has no rows: {}".format(path))
    return gene_kmers, gene_roles


def collect_sample_dbs(args):
    items = []
    if args.db:
        prefix = normalize_db_prefix(args.db)
        sample = args.sample if args.sample else infer_sample_name(prefix)
        validate_safe_component(sample, "sample")
        items.append((sample, prefix))
    elif args.db_dir:
        db_dir = Path(args.db_dir)
        if not db_dir.is_dir():
            raise FileNotFoundError("DB directory not found: {}".format(db_dir))
        pres = sorted(db_dir.glob("*.kmc_pre"))
        if not pres:
            raise ValueError("No *.kmc_pre files found in {}".format(db_dir))
        for pre in pres:
            prefix = normalize_db_prefix(pre)
            sample = infer_sample_name(prefix)
            validate_safe_component(sample, "sample")
            items.append((sample, prefix))
    elif args.db_list:
        list_path = Path(args.db_list)
        if not list_path.is_file():
            raise FileNotFoundError("DB list file not found: {}".format(list_path))
        base_dir = list_path.resolve().parent
        with list_path.open() as handle:
            first_data = True
            for line_no, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if (
                    first_data
                    and len(cols) >= 2
                    and cols[0].lower() == "sample"
                    and cols[1].lower() in ("db_prefix", "db", "kmc_db", "path")
                ):
                    first_data = False
                    continue
                first_data = False
                if len(cols) == 1:
                    raw_prefix = Path(cols[0])
                    if not raw_prefix.is_absolute():
                        raw_prefix = (base_dir / raw_prefix).resolve()
                    prefix = normalize_db_prefix(raw_prefix)
                    sample = infer_sample_name(prefix)
                elif len(cols) >= 2:
                    sample = cols[0].strip()
                    raw_prefix = Path(cols[1].strip())
                    if not raw_prefix.is_absolute():
                        raw_prefix = (base_dir / raw_prefix).resolve()
                    prefix = normalize_db_prefix(raw_prefix)
                else:
                    raise ValueError("Invalid line {} in {}".format(line_no, list_path))
                validate_safe_component(sample, "sample")
                items.append((sample, prefix))
        if not items:
            raise ValueError("No KMC databases found in {}".format(list_path))

    seen_samples = set()
    for sample, prefix in items:
        if sample in seen_samples:
            raise ValueError("Duplicate sample name: {}".format(sample))
        seen_samples.add(sample)
        check_kmc_db(prefix)
    return items


def atomic_write_lines(path, lines_iter):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".tsv", dir=str(path.parent))
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as handle:
            for line in lines_iter:
                handle.write(line)
        os.replace(str(tmp_path), str(path))
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def run_query_mmap(query_mmap_bin, db_prefix, panel_fasta, raw_out=None):
    """Run query_mmap and return (counts_by_canonical_kmer, metadata_dict)."""
    start = time.monotonic()
    cmd = [str(query_mmap_bin), str(db_prefix), str(panel_fasta)]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    counts = {}
    stored_min_count = None
    counter_bytes = None
    raw_handle = None
    raw_tmp = None

    try:
        if raw_out is not None:
            raw_out = Path(raw_out)
            raw_out.parent.mkdir(parents=True, exist_ok=True)
            fd, raw_tmp_name = tempfile.mkstemp(
                prefix=".tmp_raw_", suffix=".tsv", dir=str(raw_out.parent)
            )
            raw_tmp = Path(raw_tmp_name)
            raw_handle = os.fdopen(fd, "w")

        header_seen = False
        for line_no, raw_line in enumerate(proc.stdout, start=1):
            if raw_handle is not None:
                raw_handle.write(raw_line)
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if not header_seen:
                header_seen = True
                expected_header = (
                    "kmer_id\tsequence\tstored_count\tstored_min_count\tcounter_bytes"
                )
                if line != expected_header:
                    raise RuntimeError(
                        "Unexpected query_mmap header: {!r}".format(line)
                    )
                continue
            fields = line.split("\t")
            if len(fields) != 5:
                raise RuntimeError(
                    "Malformed query_mmap line {}: {!r}".format(line_no, line)
                )
            seq = canonical_kmer(fields[1])
            count = int(fields[2])
            min_c = int(fields[3])
            c_bytes = int(fields[4])
            if stored_min_count is None:
                stored_min_count = min_c
                counter_bytes = c_bytes
            elif stored_min_count != min_c or counter_bytes != c_bytes:
                raise RuntimeError(
                    "Inconsistent KMC metadata in query_mmap output for {}".format(
                        db_prefix
                    )
                )
            if seq in counts and counts[seq] != count:
                raise RuntimeError(
                    "Conflicting counts for canonical k-mer {}: {} vs {}".format(
                        seq, counts[seq], count
                    )
                )
            counts[seq] = count

        stderr_text = proc.stderr.read()
        retcode = proc.wait()
        if retcode != 0:
            raise RuntimeError(
                "query_mmap failed (exit {}) on {}: {}".format(
                    retcode, db_prefix, stderr_text.strip()
                )
            )
        if not counts:
            raise RuntimeError("query_mmap returned no k-mer rows for {}".format(db_prefix))

        if raw_handle is not None:
            raw_handle.close()
            raw_handle = None
            os.replace(str(raw_tmp), str(raw_out))
            raw_tmp = None
    except Exception:
        if raw_handle is not None:
            raw_handle.close()
        if raw_tmp is not None and raw_tmp.exists():
            raw_tmp.unlink()
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        raise

    elapsed = time.monotonic() - start
    max_counter = (1 << (8 * counter_bytes)) - 1
    nonzero = sum(1 for v in counts.values() if v > 0)
    zeros = len(counts) - nonzero
    saturated = sum(1 for v in counts.values() if v >= max_counter)
    meta = {
        "panel_kmers": len(counts),
        "nonzero_kmers": nonzero,
        "zero_kmers": zeros,
        "stored_min_count": stored_min_count,
        "counter_bytes": counter_bytes,
        "max_counter": max_counter,
        "saturated_kmers": saturated,
        "elapsed_seconds": elapsed,
    }
    return counts, meta


def format_realcount_lines(counts, kmers_subset=None, nonzero_only=False):
    keys = sorted(counts) if kmers_subset is None else sorted(kmers_subset)
    written = 0
    lines = []
    for kmer in keys:
        if kmer not in counts:
            raise KeyError(
                "Membership k-mer {} was not present in panel query results".format(kmer)
            )
        value = counts[kmer]
        if nonzero_only and value == 0:
            continue
        lines.append("{}\t{}\n".format(kmer, value))
        written += 1
    return lines, written


def planned_output_paths(sample, args, gene_roles):
    paths = []
    if args.membership:
        out_dir = Path(args.out_dir)
        for gene_id, role in sorted(gene_roles.items()):
            if args.split_by_role:
                paths.append(out_dir / role / gene_id / (sample + args.out_suffix))
            else:
                paths.append(out_dir / gene_id / (sample + args.out_suffix))
        if args.save_raw:
            paths.append(out_dir / "raw" / (sample + ".query_mmap.tsv"))
    elif args.out:
        paths.append(Path(args.out))
        if args.raw_out:
            paths.append(Path(args.raw_out))
    else:
        out_dir = Path(args.out_dir)
        paths.append(out_dir / (sample + args.out_suffix))
        if args.save_raw:
            paths.append(out_dir / "raw" / (sample + ".query_mmap.tsv"))
    return paths


def process_one_sample(
    sample, db_prefix, query_mmap_bin, panel_fasta, args, gene_kmers, gene_roles
):
    raw_out = None
    if args.raw_out:
        raw_out = Path(args.raw_out)
    elif args.save_raw:
        raw_out = Path(args.out_dir) / "raw" / (sample + ".query_mmap.tsv")

    counts, meta = run_query_mmap(query_mmap_bin, db_prefix, panel_fasta, raw_out=raw_out)

    if meta["stored_min_count"] > 1 and meta["zero_kmers"] > 0 and not args.nonzero_only:
        sys.stderr.write(
            "WARNING: sample {} KMC DB has stored_min_count={} ({} panel k-mers have stored_count=0, "
            "which means < {} reads rather than necessarily 0 reads)\n".format(
                sample,
                meta["stored_min_count"],
                meta["zero_kmers"],
                meta["stored_min_count"],
            )
        )
    if meta["saturated_kmers"] > 0:
        sys.stderr.write(
            "WARNING: sample {} has {} k-mers at KMC counter ceiling {} (counter_bytes={})\n".format(
                sample,
                meta["saturated_kmers"],
                meta["max_counter"],
                meta["counter_bytes"],
            )
        )

    total_written = 0
    if args.membership:
        out_dir = Path(args.out_dir)
        for gene_id, kmers_set in sorted(gene_kmers.items()):
            role = gene_roles[gene_id]
            if args.split_by_role:
                target_path = out_dir / role / gene_id / (sample + args.out_suffix)
            else:
                target_path = out_dir / gene_id / (sample + args.out_suffix)
            lines, written = format_realcount_lines(
                counts, kmers_subset=kmers_set, nonzero_only=args.nonzero_only
            )
            atomic_write_lines(target_path, lines)
            total_written += written
    else:
        target_path = (
            Path(args.out)
            if args.out
            else Path(args.out_dir) / (sample + args.out_suffix)
        )
        lines, written = format_realcount_lines(
            counts, kmers_subset=None, nonzero_only=args.nonzero_only
        )
        atomic_write_lines(target_path, lines)
        total_written = written

    return {
        "sample": sample,
        "db_prefix": str(db_prefix),
        "panel_kmers": meta["panel_kmers"],
        "written_kmers": total_written,
        "nonzero_kmers": meta["nonzero_kmers"],
        "zero_kmers": meta["zero_kmers"],
        "stored_min_count": meta["stored_min_count"],
        "counter_bytes": meta["counter_bytes"],
        "saturated_kmers": meta["saturated_kmers"],
        "elapsed_seconds": "{:.4f}".format(meta["elapsed_seconds"]),
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Query a canonical 31-mer sample reads KMC database directly against "
            "a panel FASTA/TSV using read-only mmap (`query_mmap`) and write "
            "TKCN-compatible headerless `*.realcount.tsv` count tables."
        )
    )
    db_group = parser.add_mutually_exclusive_group(required=True)
    db_group.add_argument(
        "--db",
        help="Single sample reads KMC database prefix (without .kmc_pre/.kmc_suf).",
    )
    db_group.add_argument(
        "--db-dir",
        help="Directory containing sample *.kmc_pre/*.kmc_suf databases (batch mode).",
    )
    db_group.add_argument(
        "--db-list",
        help=(
            "Text/TSV list of sample KMC DB prefixes (1 column: db_prefix; "
            "or 2 columns: sample<TAB>db_prefix)."
        ),
    )

    parser.add_argument(
        "--panel",
        required=True,
        help="Panel FASTA (e.g. TaFTB1_7B.panel.fa or all_panels.k31.fa) or panel TSV.",
    )
    parser.add_argument(
        "-o",
        "--out",
        help="Output *.realcount.tsv path (single --db mode without --membership).",
    )
    parser.add_argument(
        "--out-dir",
        help="Output directory (for batch mode, --membership splitting, or single --db).",
    )
    parser.add_argument(
        "--out-suffix",
        default=".realcount.tsv",
        help="Output filename suffix in --out-dir mode (default: .realcount.tsv).",
    )
    parser.add_argument(
        "--sample",
        help="Sample name override in single --db mode.",
    )
    parser.add_argument(
        "--membership",
        help=(
            "Optional membership.tsv (kmer, gene_id[, role]) to split a joint panel "
            "into <out-dir>/<gene_id>/<sample><out-suffix>."
        ),
    )
    parser.add_argument(
        "--split-by-role",
        action="store_true",
        help="With --membership, write under <out-dir>/<role>/<gene_id>/.",
    )
    parser.add_argument(
        "--nonzero-only",
        action="store_true",
        help=(
            "Omit k-mers with stored_count == 0 (matches legacy kmc_tools simple "
            "intersect -ocleft + dump -s). By default, all panel k-mers including "
            "0 counts are written."
        ),
    )
    parser.add_argument(
        "--raw-out",
        help="Save raw 5-column query_mmap TSV output (single --db mode).",
    )
    parser.add_argument(
        "--save-raw",
        action="store_true",
        help="Save raw 5-column query_mmap TSV files under <out-dir>/raw/.",
    )
    parser.add_argument(
        "--summary",
        help=(
            "Write per-sample query summary TSV. In batch mode (--db-dir/--db-list), "
            "defaults to <out-dir>/query_summary.tsv."
        ),
    )
    parser.add_argument(
        "--query-mmap",
        help=(
            "Path to query_mmap executable. Overrides QUERY_MMAP_BIN, "
            "src/query_mmap/query_mmap, PATH, and default PKU-IAAS path."
        ),
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=4,
        help="Parallel worker threads for batch mode (default: 4).",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.workers < 1:
        raise SystemExit("ERROR: --workers must be >= 1")
    if args.split_by_role and not args.membership:
        raise SystemExit("ERROR: --split-by-role requires --membership")

    batch_mode = bool(args.db_dir or args.db_list)
    if batch_mode:
        if args.out or args.raw_out or args.sample:
            raise SystemExit(
                "ERROR: --out, --raw-out, and --sample are only valid with single --db mode"
            )
        if not args.out_dir:
            raise SystemExit("ERROR: --out-dir is required in batch mode")
    else:
        if args.membership:
            if not args.out_dir or args.out:
                raise SystemExit(
                    "ERROR: --membership requires --out-dir (and cannot be combined with --out)"
                )
        else:
            if bool(args.out) == bool(args.out_dir):
                raise SystemExit(
                    "ERROR: specify exactly one of --out or --out-dir in single --db mode"
                )
        if args.save_raw and not args.out_dir:
            raise SystemExit("ERROR: --save-raw requires --out-dir")

    query_mmap_bin = resolve_query_mmap(args.query_mmap)
    samples = collect_sample_dbs(args)

    gene_kmers, gene_roles = ({}, {})
    if args.membership:
        gene_kmers, gene_roles = load_membership(args.membership)

    summary_path = None
    if args.summary:
        summary_path = Path(args.summary)
    elif batch_mode:
        summary_path = Path(args.out_dir) / "query_summary.tsv"

    if not args.force:
        existing = []
        for sample, _ in samples:
            for out_path in planned_output_paths(sample, args, gene_roles):
                if out_path.exists():
                    existing.append(out_path)
        if summary_path is not None and summary_path.exists():
            existing.append(summary_path)
        if existing:
            preview = ", ".join(str(p) for p in existing[:5])
            raise SystemExit(
                "ERROR: output already exists ({}{}). Use --force to overwrite.".format(
                    preview, " ..." if len(existing) > 5 else ""
                )
            )

    with tempfile.TemporaryDirectory(prefix="tkcn_query_mmap_") as tmp_dir:
        panel_fasta = prepare_panel_fasta(args.panel, tmp_dir)

        def worker(item):
            sample, db_prefix = item
            return process_one_sample(
                sample,
                db_prefix,
                query_mmap_bin,
                panel_fasta,
                args,
                gene_kmers,
                gene_roles,
            )

        results = []
        if len(samples) == 1 or args.workers == 1:
            for idx, item in enumerate(samples, start=1):
                results.append(worker(item))
                if len(samples) > 1 and (idx % 20 == 0 or idx == len(samples)):
                    sys.stderr.write(
                        "Queried {}/{} samples\n".format(idx, len(samples))
                    )
        else:
            with ThreadPoolExecutor(max_workers=min(args.workers, len(samples))) as pool:
                for idx, row in enumerate(pool.map(worker, samples), start=1):
                    results.append(row)
                    if idx % 20 == 0 or idx == len(samples):
                        sys.stderr.write(
                            "Queried {}/{} samples\n".format(idx, len(samples))
                        )

    if summary_path is not None:
        fields = [
            "sample",
            "db_prefix",
            "panel_kmers",
            "written_kmers",
            "nonzero_kmers",
            "zero_kmers",
            "stored_min_count",
            "counter_bytes",
            "saturated_kmers",
            "elapsed_seconds",
        ]
        summary_lines = ["\t".join(fields) + "\n"]
        for row in results:
            summary_lines.append(
                "\t".join(str(row[col]) for col in fields) + "\n"
            )
        atomic_write_lines(summary_path, summary_lines)
        sys.stderr.write("Summary written to: {}\n".format(summary_path))

    if len(results) == 1:
        r = results[0]
        sys.stderr.write(
            "Queried sample {} ({} panel k-mers: {} nonzero, {} zero, min_count={}, counter_bytes={}, {}s)\n".format(
                r["sample"],
                r["panel_kmers"],
                r["nonzero_kmers"],
                r["zero_kmers"],
                r["stored_min_count"],
                r["counter_bytes"],
                r["elapsed_seconds"],
            )
        )


if __name__ == "__main__":
    main()
