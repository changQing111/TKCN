#!/usr/bin/env python3
"""Filter KMC-dumped target-specific k-mers into a counting panel."""

import argparse
import math
import os
import sys


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Quality-filter target-specific k-mers dumped by KMC and write "
            "a panel TSV/FASTA."
        )
    )
    parser.add_argument("input_tsv", help="KMC dump TSV: kmer count, with optional header")
    parser.add_argument(
        "--out-prefix",
        default="target_specific_panel",
        help="Output prefix for .tsv, .fa, and .reject.tsv",
    )
    parser.add_argument(
        "--min-gc",
        type=float,
        default=35.0,
        help="Minimum GC percent to keep, default: 35.0",
    )
    parser.add_argument(
        "--max-gc",
        type=float,
        default=70.0,
        help="Maximum GC percent to keep, default: 70.0",
    )
    parser.add_argument(
        "--max-homopolymer",
        type=int,
        default=4,
        help="Maximum allowed homopolymer run length, default: 4",
    )
    parser.add_argument(
        "--min-entropy",
        type=float,
        default=1.6,
        help="Minimum Shannon entropy over A/C/G/T, default: 1.6",
    )
    parser.add_argument(
        "--max-repeat-period",
        type=int,
        default=6,
        help="Reject perfect simple repeats with motif period <= this value, default: 6",
    )
    parser.add_argument(
        "--anchor-fasta",
        help="Optional target copy FASTA used to annotate positions and thin dense k-mers",
    )
    parser.add_argument(
        "--min-spacing",
        type=int,
        default=0,
        help="Minimum distance between kept k-mer starts on the same anchor sequence",
    )
    parser.add_argument(
        "--id-prefix",
        default="panel",
        help="FASTA record prefix, default: panel",
    )
    return parser.parse_args()


def ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def iter_kmers(path):
    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("line {} has fewer than 2 columns: {}".format(line_no, line))
            if fields[0].lower() == "kmer":
                continue
            yield fields[0].upper(), fields[1]


def gc_percent(kmer):
    return (100.0 * (kmer.count("G") + kmer.count("C"))) / len(kmer)


def max_homopolymer_run(kmer):
    longest = 0
    current = 0
    previous = None
    for base in kmer:
        if base == previous:
            current += 1
        else:
            current = 1
            previous = base
        longest = max(longest, current)
    return longest


def shannon_entropy(kmer):
    entropy = 0.0
    length = float(len(kmer))
    for base in "ACGT":
        count = kmer.count(base)
        if count:
            p = count / length
            entropy -= p * math.log(p, 2)
    return entropy


def repeat_period(kmer, max_period):
    for period in range(1, min(max_period, len(kmer)) + 1):
        motif = kmer[:period]
        repeated = (motif * ((len(kmer) // period) + 1))[: len(kmer)]
        if repeated == kmer:
            return period
    return None


def reject_reason(kmer, args):
    if any(base not in "ACGT" for base in kmer):
        return "contains_non_acgt"

    gc = gc_percent(kmer)
    if gc < args.min_gc:
        return "gc_lt_{}".format(format_threshold(args.min_gc))
    if gc > args.max_gc:
        return "gc_gt_{}".format(format_threshold(args.max_gc))

    max_run = max_homopolymer_run(kmer)
    if max_run > args.max_homopolymer:
        return "homopolymer_gt_{}".format(args.max_homopolymer)

    period = repeat_period(kmer, args.max_repeat_period)
    if period:
        return "simple_repeat_period_{}".format(period)

    entropy = shannon_entropy(kmer)
    if entropy < args.min_entropy:
        return "entropy_lt_{}".format(format_threshold(args.min_entropy))

    return None


def format_threshold(value):
    if int(value) == value:
        return str(int(value))
    return str(value)


def read_fasta(path):
    records = []
    name = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    return records


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def find_anchor_position(kmer, records):
    rc = reverse_complement(kmer)
    for seq_id, seq in records:
        pos = seq.find(kmer)
        if pos != -1:
            return seq_id, pos + 1, pos + len(kmer), "+"
        if rc != kmer:
            pos = seq.find(rc)
            if pos != -1:
                return seq_id, pos + 1, pos + len(kmer), "-"
    return ".", ".", ".", "."


def apply_spacing(records, min_spacing):
    if min_spacing <= 0:
        return records, []

    kept = []
    rejected = []
    last_start_by_anchor = {}
    for record in sorted(records, key=spacing_sort_key):
        anchor_id = record["anchor_id"]
        anchor_start = record["anchor_start"]
        if anchor_id == "." or anchor_start == ".":
            kept.append(record)
            continue

        last_start = last_start_by_anchor.get(anchor_id)
        if last_start is not None and anchor_start - last_start < min_spacing:
            rejected.append((record["kmer"], record["count"], "spacing_lt_{}".format(min_spacing)))
            continue

        kept.append(record)
        last_start_by_anchor[anchor_id] = anchor_start

    return kept, rejected


def spacing_sort_key(record):
    anchor_id = record["anchor_id"]
    anchor_start = record["anchor_start"]
    if anchor_id == "." or anchor_start == ".":
        return ("~", sys.maxsize, record["kmer"])
    return (anchor_id, anchor_start, record["kmer"])


def build_panel(args):
    anchor_records = read_fasta(args.anchor_fasta) if args.anchor_fasta else []
    accepted = []
    rejected = []

    for kmer, count in iter_kmers(args.input_tsv):
        reason = reject_reason(kmer, args)
        if reason:
            rejected.append((kmer, count, reason))
            continue

        anchor_id, anchor_start, anchor_end, strand = find_anchor_position(kmer, anchor_records)
        accepted.append(
            {
                "kmer": kmer,
                "count": count,
                "gc_percent": gc_percent(kmer),
                "max_homopolymer": max_homopolymer_run(kmer),
                "entropy": shannon_entropy(kmer),
                "anchor_id": anchor_id,
                "anchor_start": anchor_start,
                "anchor_end": anchor_end,
                "strand": strand,
            }
        )

    accepted, spacing_rejected = apply_spacing(accepted, args.min_spacing)
    rejected.extend(spacing_rejected)
    return accepted, rejected


def write_outputs(args, accepted, rejected):
    panel_tsv = args.out_prefix + ".tsv"
    panel_fa = args.out_prefix + ".fa"
    reject_tsv = args.out_prefix + ".reject.tsv"

    for path in (panel_tsv, panel_fa, reject_tsv):
        ensure_parent_dir(path)

    with open(panel_tsv, "w") as out:
        out.write(
            "\t".join(
                [
                    "kmer",
                    "count",
                    "gc_percent",
                    "max_homopolymer",
                    "entropy",
                    "anchor_id",
                    "anchor_start",
                    "anchor_end",
                    "strand",
                ]
            )
            + "\n"
        )
        for record in accepted:
            out.write(
                "{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{}\t{}\n".format(
                    record["kmer"],
                    record["count"],
                    record["gc_percent"],
                    record["max_homopolymer"],
                    record["entropy"],
                    record["anchor_id"],
                    record["anchor_start"],
                    record["anchor_end"],
                    record["strand"],
                )
            )

    with open(panel_fa, "w") as out:
        for index, record in enumerate(accepted, start=1):
            out.write(
                ">{}_{:06d} count={} gc={:.2f} max_homopolymer={}\n{}\n".format(
                    args.id_prefix,
                    index,
                    record["count"],
                    record["gc_percent"],
                    record["max_homopolymer"],
                    record["kmer"],
                )
            )

    with open(reject_tsv, "w") as out:
        out.write("kmer\tcount\treason\n")
        for kmer, count, reason in rejected:
            out.write("{}\t{}\t{}\n".format(kmer, count, reason))

    return panel_tsv, panel_fa, reject_tsv


def main():
    args = parse_args()
    accepted, rejected = build_panel(args)
    panel_tsv, panel_fa, reject_tsv = write_outputs(args, accepted, rejected)
    sys.stderr.write("Kept {} k-mers; rejected {} k-mers\n".format(len(accepted), len(rejected)))
    sys.stderr.write("Panel TSV: {}\n".format(panel_tsv))
    sys.stderr.write("Panel FASTA: {}\n".format(panel_fa))
    sys.stderr.write("Reject TSV: {}\n".format(reject_tsv))


if __name__ == "__main__":
    main()
