"""Run: python3 -m unittest discover -s test -p 'test_query_panel_kmc_mmap.py' -v"""
import csv
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "src" / "query_panel_kmc_mmap.py"
spec = importlib.util.spec_from_file_location("query_panel_kmc_mmap", str(SCRIPT))
qmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qmod)


def find_bin(env_name, fallback_path, cmd_name):
    candidates = [os.environ.get(env_name), fallback_path, shutil.which(cmd_name)]
    for c in candidates:
        if c and Path(c).is_file() and os.access(str(c), os.X_OK):
            return str(c)
    return None


KMC_BIN = find_bin("KMC_BIN", "/data/changq/bio_soft/KMC3.2.4/bin/kmc", "kmc")
KMC_TOOLS_BIN = find_bin(
    "KMC_TOOLS_BIN", "/data/changq/bio_soft/KMC3.2.4/bin/kmc_tools", "kmc_tools"
)
PYTHON_SCIPY = find_bin("PYTHON", "/data/changq/miniforge3/bin/python", "python3")

try:
    QUERY_MMAP_BIN = str(qmod.resolve_query_mmap())
except FileNotFoundError:
    QUERY_MMAP_BIN = None


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_unique_31mers(n, seed=42):
    # Deterministic distinct 31-mers with no collisions across canonical forms.
    bases = "ACGT"
    out = []
    seen = set()
    x = seed
    while len(out) < n:
        chars = []
        for _ in range(31):
            x = (1103515245 * x + 12345) & 0x7FFFFFFF
            chars.append(bases[(x >> 16) & 3])
        seq = "".join(chars)
        can = qmod.canonical_kmer(seq)
        if can not in seen:
            seen.add(can)
            out.append(seq)
    return out


class QueryPanelUnitTests(unittest.TestCase):
    def test_canonical_and_name_helpers(self):
        seq = "TTTTACGTCAGTACGATCGTACCGATGCTAG"
        can = qmod.canonical_kmer(seq)
        self.assertEqual(can, min(seq, rc(seq)))
        self.assertEqual(qmod.canonical_kmer(rc(seq).lower()), can)
        self.assertEqual(qmod.infer_sample_name("/tmp/S01.reads.k31.kmc_pre"), "S01")
        self.assertEqual(qmod.infer_sample_name("/tmp/S02.k31"), "S02")
        with self.assertRaises(ValueError):
            qmod.canonical_kmer("ACGT")
        with self.assertRaises(ValueError):
            qmod.validate_safe_component("../bad", "sample")

    def test_prepare_panel_fasta_from_tsv_and_fasta(self):
        kmers = make_unique_31mers(4, seed=7)
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fa = root / "panel.fa"
            fa.write_text(
                "".join(">k{}\n{}\n".format(i, s) for i, s in enumerate(kmers, 1))
            )
            self.assertEqual(qmod.prepare_panel_fasta(fa, d), fa)

            headered_tsv = root / "panel.tsv"
            headered_tsv.write_text(
                "kmer\tcount\tgc\n"
                + "".join("{}\t1\t50.0\n".format(s) for s in kmers)
                + "{}\t1\t50.0\n".format(rc(kmers[0]))
            )
            converted = qmod.prepare_panel_fasta(headered_tsv, d)
            lines = [l for l in converted.read_text().splitlines() if not l.startswith(">")]
            self.assertEqual(lines, sorted(qmod.canonical_kmer(s) for s in kmers))

            headerless_tsv = root / "shared.tsv"
            headerless_tsv.write_text("".join("{}\t1\n".format(s) for s in kmers[:2]))
            converted2 = qmod.prepare_panel_fasta(headerless_tsv, d)
            lines2 = [l for l in converted2.read_text().splitlines() if not l.startswith(">")]
            self.assertEqual(lines2, sorted(qmod.canonical_kmer(s) for s in kmers[:2]))

    def test_load_membership_and_role_conflict(self):
        kmers = make_unique_31mers(3, seed=11)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "membership.tsv"
            p.write_text(
                "panel_id\tkmer\tgene_id\trole\n"
                "p1\t{0}\tG1\tsingle_copy\n"
                "p2\t{1}\tT1\ttarget\n"
                "p2\t{2}\tT2\ttarget\n".format(kmers[0], kmers[1], rc(kmers[1]))
            )
            gk, gr = qmod.load_membership(p)
            self.assertEqual(gr, {"G1": "single_copy", "T1": "target", "T2": "target"})
            self.assertEqual(gk["T1"], gk["T2"])


@unittest.skipUnless(
    KMC_BIN and KMC_TOOLS_BIN and QUERY_MMAP_BIN,
    "Requires kmc, kmc_tools, and query_mmap executables",
)
class QueryPanelIntegrationTests(unittest.TestCase):
    def _build_kmc_db(self, reads_fa, db_prefix, tmp_dir, ci=1, cs=1000000000):
        kmc_tmp = Path(tmp_dir) / ("kmc_tmp_" + Path(db_prefix).name)
        kmc_tmp.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                KMC_BIN,
                "-k31",
                "-ci{}".format(ci),
                "-cs{}".format(cs),
                "-fm",
                "-t2",
                str(reads_fa),
                str(db_prefix),
                str(kmc_tmp),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def test_exact_parity_with_kmc_intersect_and_downstream_nb(self):
        all_kmers = make_unique_31mers(24, seed=101)
        single_a = all_kmers[0:6]
        single_b = all_kmers[6:12]
        target = all_kmers[12:24]

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            # Sample S2: single-copy depth ~10x, target copy number = 2 (~20x),
            # with one absent target k-mer (0) and one singleton (1) for edge testing.
            s2_counts = {}
            for i, s in enumerate(single_a):
                s2_counts[s] = 9 + (i % 3)  # 9,10,11,9,10,11 (mean 10)
            for i, s in enumerate(single_b):
                s2_counts[s] = 10 + (i % 2)  # 10,11,...
            target_s2_vals = [20, 21, 19, 20, 20, 21, 19, 20, 20, 19, 1, 0]
            for s, v in zip(target, target_s2_vals):
                s2_counts[s] = v

            # Sample S0: single-copy depth ~10x, target copy number = 0 (all 0).
            s0_counts = {}
            for i, s in enumerate(single_a + single_b):
                s0_counts[s] = 10
            for s in target:
                s0_counts[s] = 0

            def write_reads(path, count_map):
                idx = 0
                with path.open("w") as h:
                    for seq, cnt in count_map.items():
                        for rep in range(cnt):
                            idx += 1
                            # Alternate forward and reverse-complement reads to verify canonical matching.
                            out_seq = rc(seq) if (rep % 2 == 1) else seq
                            h.write(">r{}\n{}\n".format(idx, out_seq))

            db_dir = root / "dbs"
            db_dir.mkdir()
            write_reads(root / "S2.reads.fa", s2_counts)
            write_reads(root / "S0.reads.fa", s0_counts)
            self._build_kmc_db(root / "S2.reads.fa", db_dir / "S2.reads.k31", root)
            self._build_kmc_db(root / "S0.reads.fa", db_dir / "S0.reads.k31", root)

            # Build target panel FASTA (mix canonical and non-canonical orientations).
            target_panel_fa = root / "TaFTB1_7B.panel.fa"
            with target_panel_fa.open("w") as h:
                for i, seq in enumerate(target, start=1):
                    h.write(
                        ">TaFTB1_7B_{:06d} count=1 gc=50.00\n{}\n".format(
                            i, rc(seq) if i % 2 == 0 else seq
                        )
                    )

            # 1) Legacy kmc_tools simple intersect -ocleft + dump -s
            panel_db = root / "TaFTB1_7B.panel.k31"
            self._build_kmc_db(target_panel_fa, panel_db, root, ci=1, cs=255)
            inter_db = root / "S2.target.k31"
            legacy_tsv = root / "S2.legacy.realcount.tsv"
            subprocess.run(
                [
                    KMC_TOOLS_BIN,
                    "simple",
                    str(db_dir / "S2.reads.k31"),
                    str(panel_db),
                    "intersect",
                    str(inter_db),
                    "-ocleft",
                    "-cs1000000000",
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            subprocess.run(
                [KMC_TOOLS_BIN, "transform", str(inter_db), "dump", "-s", str(legacy_tsv)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # 2) New query_panel_kmc_mmap.py (--nonzero-only) must match legacy byte-for-byte.
            mmap_nonzero = root / "S2.mmap_nonzero.realcount.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--db",
                    str(db_dir / "S2.reads.k31"),
                    "--panel",
                    str(target_panel_fa),
                    "--out",
                    str(mmap_nonzero),
                    "--nonzero-only",
                    "--query-mmap",
                    QUERY_MMAP_BIN,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(mmap_nonzero.read_text(), legacy_tsv.read_text())

            # 3) Default mode keeps zeros + writes raw 5-col output; refuse overwrite without --force.
            mmap_full = root / "S2.mmap_full.realcount.tsv"
            raw_tsv = root / "S2.raw.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--db",
                    str(db_dir / "S2.reads.k31"),
                    "--panel",
                    str(target_panel_fa),
                    "--out",
                    str(mmap_full),
                    "--raw-out",
                    str(raw_tsv),
                    "--query-mmap",
                    QUERY_MMAP_BIN,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            full_rows = [l.split("\t") for l in mmap_full.read_text().splitlines()]
            self.assertEqual(len(full_rows), len(target))
            expected_full = sorted(
                (qmod.canonical_kmer(s), str(v)) for s, v in zip(target, target_s2_vals)
            )
            self.assertEqual([(r[0], r[1]) for r in full_rows], expected_full)
            self.assertTrue(raw_tsv.read_text().startswith("kmer_id\tsequence\tstored_count\t"))

            # Verify parity with pre-existing TaFTL query_mmap binary if present on host.
            if qmod.FALLBACK_QUERY_MMAP.is_file():
                mmap_fallback = root / "S2.fallback.realcount.tsv"
                subprocess.run(
                    [
                        sys.executable,
                        str(SCRIPT),
                        "--db",
                        str(db_dir / "S2.reads.k31"),
                        "--panel",
                        str(target_panel_fa),
                        "--out",
                        str(mmap_fallback),
                        "--query-mmap",
                        str(qmod.FALLBACK_QUERY_MMAP),
                    ],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                self.assertEqual(mmap_fallback.read_text(), mmap_full.read_text())

            # Overwrite protection without --force
            fail = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--db",
                    str(db_dir / "S2.reads.k31"),
                    "--panel",
                    str(target_panel_fa),
                    "--out",
                    str(mmap_full),
                    "--query-mmap",
                    QUERY_MMAP_BIN,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertNotEqual(fail.returncode, 0)

            # 4) Batch mode + joint membership split + downstream NB copy-number estimation!
            joint_fa = root / "joint.k31.fa"
            membership_tsv = root / "membership.tsv"
            with joint_fa.open("w") as fh, membership_tsv.open("w") as mh:
                mh.write("panel_id\tkmer\tgene_id\trole\n")
                idx = 0
                for gene_id, role, seqs in [
                    ("MarkerA", "single_copy", single_a),
                    ("MarkerB", "single_copy", single_b),
                    ("TaFTB1_7B", "target", target),
                ]:
                    for seq in seqs:
                        idx += 1
                        pid = "joint_{:06d}".format(idx)
                        can = qmod.canonical_kmer(seq)
                        fh.write(">{}\n{}\n".format(pid, can))
                        mh.write("{}\t{}\t{}\t{}\n".format(pid, can, gene_id, role))

            split_out = root / "split_counts"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--db-dir",
                    str(db_dir),
                    "--panel",
                    str(joint_fa),
                    "--membership",
                    str(membership_tsv),
                    "--split-by-role",
                    "--out-dir",
                    str(split_out),
                    "--workers",
                    "2",
                    "--query-mmap",
                    QUERY_MMAP_BIN,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertTrue((split_out / "query_summary.tsv").is_file())
            self.assertTrue((split_out / "single_copy" / "MarkerA" / "S2.realcount.tsv").is_file())
            self.assertTrue((split_out / "target" / "TaFTB1_7B" / "S0.realcount.tsv").is_file())

            # Also verify compatibility with downstream NB scripts if SciPy python is present.
            if PYTHON_SCIPY:
                check_scipy = subprocess.run(
                    [PYTHON_SCIPY, "-c", "import scipy"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                if check_scipy.returncode == 0:
                    nb_script = REPO_ROOT / "src" / "estimate_depth_from_realcount_nb.py"
                    copy_script = REPO_ROOT / "src" / "estimate_copy_from_realcount_depth.sh"
                    marker_files = sorted(
                        str(p) for p in (split_out / "single_copy").glob("*/*.realcount.tsv")
                    )
                    marker_nb = root / "single_copy_marker_nb.tsv"
                    joint_depth = root / "single_copy_depth_nb_joint.tsv"
                    copy_est = root / "TaFTB1_7B_copy_estimate.tsv"

                    r1 = subprocess.run(
                        [PYTHON_SCIPY, str(nb_script), "marker"] + marker_files,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True,
                    )
                    marker_nb.write_text(r1.stdout)
                    r2 = subprocess.run(
                        [PYTHON_SCIPY, str(nb_script), "joint", "--marker-stats", str(marker_nb)],
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True,
                    )
                    joint_depth.write_text(r2.stdout)

                    env = os.environ.copy()
                    env["PYTHON"] = PYTHON_SCIPY
                    subprocess.run(
                        [
                            "bash",
                            str(copy_script),
                            "--target",
                            str(split_out / "target" / "TaFTB1_7B"),
                            "--depth",
                            str(joint_depth),
                            "--target-method",
                            "nb",
                            "--copy-call-method",
                            "nb_likelihood",
                            "--out",
                            str(copy_est),
                        ],
                        check=True,
                        env=env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    with copy_est.open() as fh:
                        calls = {
                            row["sample"]: row["copy_call"]
                            for row in csv.DictReader(fh, delimiter="\t")
                        }
                    self.assertEqual(calls["S0"], "0")
                    self.assertEqual(calls["S2"], "2")


if __name__ == "__main__":
    unittest.main()
