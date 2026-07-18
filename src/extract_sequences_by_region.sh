#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  extract_sequences_by_region.sh -f reference.fa -r regions.txt -o output.fa [-t threads]

Description:
  Extract sequences from a FASTA file using a region file with one region per line:
  chromosome:start-end

Options:
  -f FILE   Reference FASTA file.
  -r FILE   Region file, using 1-based inclusive coordinates.
  -o FILE   Output FASTA file. The filename can be chosen freely.
  -t INT    Number of seqkit threads (default: 4).
  -h        Show this help message.

Example:
  extract_sequences_by_region.sh \
    -f reference.fa \
    -r kmer_intersect/TraesCS5A03G0935400/S4185_region.txt \
    -o S4185_TraesCS5A03G0935400.fa

  # Set SEQKIT_BIN when seqkit is not in PATH:
  SEQKIT_BIN=/path/to/seqkit extract_sequences_by_region.sh \
    -f reference.fa -r regions.txt -o extracted.fa
EOF
}

reference_fasta=""
region_file=""
output_file=""
threads=4

while getopts ":f:r:o:t:h" opt; do
    case "$opt" in
        f) reference_fasta="$OPTARG" ;;
        r) region_file="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument." >&2
            usage >&2
            exit 2
            ;;
        \?)
            echo "Error: unknown option -$OPTARG." >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$reference_fasta" || -z "$region_file" || -z "$output_file" ]]; then
    echo "Error: -f, -r and -o are required." >&2
    usage >&2
    exit 2
fi

if [[ ! -f "$reference_fasta" ]]; then
    echo "Error: reference FASTA does not exist: $reference_fasta" >&2
    exit 1
fi

if [[ ! -f "$region_file" ]]; then
    echo "Error: region file does not exist: $region_file" >&2
    exit 1
fi

if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer: $threads" >&2
    exit 2
fi

seqkit_bin="${SEQKIT_BIN:-seqkit}"
if ! command -v "$seqkit_bin" >/dev/null 2>&1; then
    echo "Error: seqkit was not found: $seqkit_bin" >&2
    echo "Set SEQKIT_BIN to the seqkit executable path." >&2
    exit 127
fi

output_dir="$(dirname -- "$output_file")"
if [[ ! -d "$output_dir" ]]; then
    echo "Error: output directory does not exist: $output_dir" >&2
    exit 1
fi

bed_file="$(mktemp)"
trap 'rm -f "$bed_file"' EXIT

# Convert 1-based inclusive regions to 0-based half-open BED intervals.
awk '
BEGIN {
    OFS = "\t"
    region_count = 0
    has_error = 0
}
{
    line = $0
    gsub(/\r$/, "", line)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
    if (line == "") {
        next
    }

    field_count = split(line, fields, /:|-/)
    for (i = 1; i <= field_count; i++) {
        gsub(/[[:space:]]/, "", fields[i])
    }

    if (field_count != 3 ||
        fields[1] == "" ||
        fields[2] !~ /^[0-9]+$/ ||
        fields[3] !~ /^[0-9]+$/ ||
        fields[2] < 1 ||
        fields[3] < fields[2]) {
        printf "Invalid region on line %d: %s\n", NR, line > "/dev/stderr"
        has_error = 1
        next
    }

    print fields[1], fields[2] - 1, fields[3], "region_" ++region_count
}
END {
    if (has_error || region_count == 0) {
        exit 1
    }
}
' "$region_file" > "$bed_file"

"$seqkit_bin" subseq \
    --bed "$bed_file" \
    --threads "$threads" \
    --update-faidx \
    --out-file "$output_file" \
    "$reference_fasta"

echo "Extracted sequences written to: $output_file"
