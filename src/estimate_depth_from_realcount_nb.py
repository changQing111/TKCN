#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

try:
    from scipy.special import gammaln
    from scipy.optimize import minimize_scalar
except ImportError as exc:
    raise SystemExit(f"ERROR: SciPy is required for NB calculations: {exc}")


class EmptyKmerCountFile(Exception):
    pass


def format_number(value):
    if value is None:
        return "NA"
    return f"{value:.10g}"


def sample_from_path(path):
    sample = Path(path).name
    for suffix in (".tsv", "_inter.k31", ".inter.k31", ".realcount", ".panel"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
    return sample


def read_counts(path):
    counts = []
    with open(path, newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise SystemExit(
                    f"ERROR: expected at least 2 columns in {path} at line {line_number}"
                )
            try:
                counts.append(float(fields[1]))
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: non-numeric count in {path} at line {line_number}"
                ) from exc
    if not counts:
        raise EmptyKmerCountFile(path)
    return counts


def poisson_loglik(counts, mu):
    if mu <= 0:
        return None
    return sum(count * math.log(mu) - mu - gammaln(count + 1) for count in counts)


def poisson_loglik_constant(counts):
    return sum(-gammaln(count + 1) for count in counts)


def nb_loglik(counts, mu, size):
    if mu <= 0 or size <= 0:
        return None
    return sum(
        gammaln(count + size)
        - gammaln(size)
        - gammaln(count + 1)
        + size * math.log(size / (size + mu))
        + count * math.log(mu / (size + mu))
        for count in counts
    )


def nb_loglik_constant(counts, size):
    if size <= 0:
        return None
    return sum(
        gammaln(count + size) - gammaln(size) - gammaln(count + 1)
        for count in counts
    )


def marker_summary(path, copy_number):
    counts = read_counts(path)
    kmer_count = len(counts)
    sum_count = sum(counts)
    mean_count = sum_count / kmer_count
    nb_depth_estimate = mean_count / copy_number
    row = {
        "sample": sample_from_path(path),
        "input_file": str(path),
        "kmer_count": kmer_count,
        "copy_number": copy_number,
        "sum_count": sum_count,
        "mean_count": mean_count,
        "var_count": None,
        "nb_depth_estimate": nb_depth_estimate,
        "nb_size": None,
        "nb_loglik": None,
        "nb_status": "insufficient_kmers",
        "nb_loglik_constant": None,
    }
    if kmer_count < 2:
        return row

    var_count = statistics.variance(counts)
    row["var_count"] = var_count
    if var_count <= mean_count:
        row["nb_status"] = "poisson_limit"
        row["nb_loglik_constant"] = poisson_loglik_constant(counts)
        row["nb_loglik"] = poisson_loglik(counts, mean_count)
        return row

    size = (mean_count * mean_count) / (var_count - mean_count)
    row["nb_status"] = "ok"
    row["nb_size"] = size
    row["nb_loglik_constant"] = nb_loglik_constant(counts, size)
    row["nb_loglik"] = nb_loglik(counts, mean_count, size)
    return row


def write_marker_rows(paths, copy_number):
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "sample",
            "input_file",
            "kmer_count",
            "copy_number",
            "sum_count",
            "mean_count",
            "var_count",
            "nb_depth_estimate",
            "nb_size",
            "nb_loglik",
            "nb_status",
            "nb_loglik_constant",
        ]
    )
    for path in paths:
        try:
            row = marker_summary(path, copy_number)
        except EmptyKmerCountFile:
            print(f"WARNING: skipped empty marker file: {path}", file=sys.stderr)
            continue
        writer.writerow(
            [
                row["sample"],
                row["input_file"],
                row["kmer_count"],
                format_number(row["copy_number"]),
                format_number(row["sum_count"]),
                format_number(row["mean_count"]),
                format_number(row["var_count"]),
                format_number(row["nb_depth_estimate"]),
                format_number(row["nb_size"]),
                format_number(row["nb_loglik"]),
                row["nb_status"],
                format_number(row["nb_loglik_constant"]),
            ]
        )


def read_marker_stats(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"ERROR: no marker NB rows in {path}")
    return rows


def marker_is_usable(row):
    return row["nb_status"] in {"ok", "poisson_limit"}


def joint_loglik(depth, markers, copy_number):
    if depth <= 0:
        return float("-inf")

    mu = copy_number * depth
    total = 0.0
    for row in markers:
        kmer_count = int(row["kmer_count"])
        sum_count = float(row["sum_count"])
        constant = float(row["nb_loglik_constant"])
        if row["nb_status"] == "ok":
            size = float(row["nb_size"])
            total += (
                constant
                + sum_count * math.log(mu / (size + mu))
                + kmer_count * size * math.log(size / (size + mu))
            )
        else:
            total += constant + sum_count * math.log(mu) - kmer_count * mu
    return total


def joint_summary(sample, rows, copy_number):
    marker_sizes = [row["nb_size"] for row in rows]
    marker_logliks = [row["nb_loglik"] for row in rows]
    marker_statuses = [row["nb_status"] for row in rows]
    usable = [row for row in rows if marker_is_usable(row)]
    usable_depths = [float(row["nb_depth_estimate"]) for row in usable]
    dropped = len(rows) - len(usable)

    result = {
        "sample": sample,
        "marker_nb_sizes": ",".join(marker_sizes),
        "marker_nb_logliks": ",".join(marker_logliks),
        "marker_nb_statuses": ",".join(marker_statuses),
        "nb_used_marker_count": str(len(usable)),
        "nb_depth_estimate": "NA",
        "nb_loglik": "NA",
        "nb_joint_status": "all_markers_failed",
        "nb_combine_method": "marker_nb_joint_mle",
    }
    if not usable:
        return result

    if all(row["nb_status"] == "poisson_limit" for row in usable):
        total_count = sum(float(row["sum_count"]) for row in usable)
        total_kmers = sum(int(row["kmer_count"]) for row in usable)
        best_depth = total_count / (total_kmers * copy_number)
    else:
        lower = 1e-8
        upper = max(1.0, 10.0 * max(usable_depths))
        optimum = minimize_scalar(
            lambda depth: -joint_loglik(depth, usable, copy_number),
            bounds=(lower, upper),
            method="bounded",
        )
        if not optimum.success or not math.isfinite(optimum.x):
            result["nb_joint_status"] = "fit_failed"
            return result
        best_depth = optimum.x

    best_loglik = joint_loglik(best_depth, usable, copy_number)
    if not math.isfinite(best_loglik):
        result["nb_joint_status"] = "fit_failed"
        return result

    result["nb_depth_estimate"] = format_number(best_depth)
    result["nb_loglik"] = format_number(best_loglik)
    result["nb_joint_status"] = "partial_marker_drop" if dropped else "ok"
    return result


def write_joint_rows(marker_stats_path, copy_number):
    rows = read_marker_stats(marker_stats_path)
    grouped = {}
    sample_order = []
    for row in rows:
        sample = row["sample"]
        if sample not in grouped:
            grouped[sample] = []
            sample_order.append(sample)
        grouped[sample].append(row)

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "sample",
            "marker_nb_sizes",
            "marker_nb_logliks",
            "marker_nb_statuses",
            "nb_used_marker_count",
            "nb_depth_estimate",
            "nb_loglik",
            "nb_joint_status",
            "nb_combine_method",
        ]
    )
    for sample in sample_order:
        row = joint_summary(sample, grouped[sample], copy_number)
        writer.writerow(
            [
                row["sample"],
                row["marker_nb_sizes"],
                row["marker_nb_logliks"],
                row["marker_nb_statuses"],
                row["nb_used_marker_count"],
                row["nb_depth_estimate"],
                row["nb_loglik"],
                row["nb_joint_status"],
                row["nb_combine_method"],
            ]
        )


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    def positive_finite(value):
        try:
            number = float(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError("must be a number") from exc
        if not math.isfinite(number) or number <= 0:
            raise argparse.ArgumentTypeError("must be finite and greater than zero")
        return number

    marker = subparsers.add_parser("marker")
    marker.add_argument("--copy-number", type=positive_finite, default=1.0)
    marker.add_argument("paths", nargs="+")

    joint = subparsers.add_parser("joint")
    joint.add_argument("--copy-number", type=positive_finite, default=1.0)
    joint.add_argument("--marker-stats", required=True)

    return parser


def main():
    args = build_parser().parse_args()
    if args.command == "marker":
        write_marker_rows(args.paths, args.copy_number)
        return 0
    if args.command == "joint":
        write_joint_rows(args.marker_stats, args.copy_number)
        return 0
    raise SystemExit(f"ERROR: unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
