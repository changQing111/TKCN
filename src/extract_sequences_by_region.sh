#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  extract_sequences_by_region.sh -f reference.fa -r regions.txt -o output.fa

Description:
  Extract sequences from a FASTA file using a region file with one region per line:
  chromosome:start-end

Options:
  -f FILE   Reference FASTA file.
  -r FILE   Region file, using 1-based inclusive coordinates.
  -o FILE   Output FASTA file. The filename can be chosen freely.
  -h        Show this help message.

Example:
  extract_sequences_by_region.sh \
    -f reference.fa \
    -r kmer_intersect/TraesCS5A03G0935400/S4185_region.txt \
    -o S4185_TraesCS5A03G0935400.fa

  # Set SAMTOOLS_BIN when samtools is not in PATH:
  SAMTOOLS_BIN=/path/to/samtools extract_sequences_by_region.sh \
    -f reference.fa -r regions.txt -o extracted.fa
EOF
}

reference_fasta=""
region_file=""
output_file=""

while getopts ":f:r:o:h" opt; do
    case "$opt" in
        f) reference_fasta="$OPTARG" ;;
        r) region_file="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
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

samtools_bin="${SAMTOOLS_BIN:-samtools}"
if ! command -v "$samtools_bin" >/dev/null 2>&1; then
    echo "Error: samtools was not found: $samtools_bin" >&2
    echo "Set SAMTOOLS_BIN to the samtools executable path." >&2
    exit 127
fi

output_dir="$(dirname -- "$output_file")"
if [[ ! -d "$output_dir" ]]; then
    echo "Error: output directory does not exist: $output_dir" >&2
    exit 1
fi

"$samtools_bin" faidx \
    -r "$region_file" \
    -o "$output_file" \
    "$reference_fasta"

echo "Extracted sequences written to: $output_file"
