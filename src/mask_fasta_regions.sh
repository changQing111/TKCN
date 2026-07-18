#!/usr/bin/env bash
set -euo pipefail

BEDTOOLS=${BEDTOOLS:-}
PYTHON=${PYTHON:-/data/changq/miniforge3/bin/python}
MASK_CHAR=${MASK_CHAR:-N}

usage() {
    cat >&2 <<'EOF'
Usage:
  mask_fasta_regions.sh -g GENOME.fa -r REGIONS.txt -o OUT.masked.fa
  mask_fasta_regions.sh GENOME.fa REGIONS.txt OUT.masked.fa

Example:
  bash mask_fasta_regions.sh \
    -g /home/changq/my_data/Project/CSv2.1_sentieon_bwa_index/iwgsc_refseqv2.1_assembly.fa \
    -r CS_all_TaFT-B1_region.txt \
    -o CS.TaFT-B1_regions.masked.fa

REGIONS.txt format:
  Chr7B:10013330-10014587
  Chr7B:10025677-10026934
  Chr7B:10038024-10039281

Input coordinates are 1-based closed intervals. The script converts them to
0-based half-open BED intervals before masking.

Environment overrides:
  BEDTOOLS=/path/to/bedtools PYTHON=/path/to/python MASK_CHAR=N bash mask_fasta_regions.sh ...
EOF
}

find_bedtools() {
    if [[ -n "$BEDTOOLS" ]]; then
        [[ -x "$BEDTOOLS" ]] || { echo "ERROR: BEDTOOLS is not executable: $BEDTOOLS" >&2; exit 1; }
        return
    fi

    local candidate
    for candidate in \
        /data/changq/miniforge3/bin/bedtools \
        /home/changq/my_data/bio_soft/cactus-bin-v2.9.8/bin/bedtools \
        /home/changq/my_data/bio_soft/cactus-bin-v2.5.2/bin/bedtools
    do
        if [[ -x "$candidate" ]]; then
            BEDTOOLS=$candidate
            return
        fi
    done

    if command -v bedtools >/dev/null 2>&1; then
        BEDTOOLS=$(command -v bedtools)
        return
    fi

    echo "ERROR: bedtools not found. Set BEDTOOLS=/path/to/bedtools" >&2
    exit 1
}

find_python() {
    if [[ -x "$PYTHON" ]]; then
        return
    fi
    if command -v python3 >/dev/null 2>&1; then
        PYTHON=$(command -v python3)
        return
    fi
    if command -v python >/dev/null 2>&1; then
        PYTHON=$(command -v python)
        return
    fi

    echo "ERROR: python not found. Set PYTHON=/path/to/python" >&2
    exit 1
}

GENOME=""
REGIONS=""
OUT_FASTA=""
KEEP_BED=""
POSITIONAL=()

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -g|--genome)
            GENOME=${2:?ERROR: --genome needs a FASTA path}
            shift 2
            ;;
        --genome=*)
            GENOME=${1#*=}
            shift
            ;;
        -r|--regions)
            REGIONS=${2:?ERROR: --regions needs a region file}
            shift 2
            ;;
        --regions=*)
            REGIONS=${1#*=}
            shift
            ;;
        -o|--output)
            OUT_FASTA=${2:?ERROR: --output needs an output FASTA path}
            shift 2
            ;;
        --output=*)
            OUT_FASTA=${1#*=}
            shift
            ;;
        --keep-bed)
            KEEP_BED=${2:?ERROR: --keep-bed needs an output BED path}
            shift 2
            ;;
        --keep-bed=*)
            KEEP_BED=${1#*=}
            shift
            ;;
        -*)
            echo "ERROR: unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

if [[ "${#POSITIONAL[@]}" -gt 3 ]]; then
    usage
    exit 1
fi
if [[ "${#POSITIONAL[@]}" -ge 1 ]]; then
    [[ -z "$GENOME" ]] || { echo "ERROR: genome given twice" >&2; exit 1; }
    GENOME=${POSITIONAL[0]}
fi
if [[ "${#POSITIONAL[@]}" -ge 2 ]]; then
    [[ -z "$REGIONS" ]] || { echo "ERROR: regions given twice" >&2; exit 1; }
    REGIONS=${POSITIONAL[1]}
fi
if [[ "${#POSITIONAL[@]}" -ge 3 ]]; then
    [[ -z "$OUT_FASTA" ]] || { echo "ERROR: output given twice" >&2; exit 1; }
    OUT_FASTA=${POSITIONAL[2]}
fi

if [[ -z "$GENOME" || -z "$REGIONS" || -z "$OUT_FASTA" ]]; then
    usage
    exit 1
fi
if [[ ! -s "$GENOME" ]]; then
    echo "ERROR: genome FASTA not found or empty: $GENOME" >&2
    exit 1
fi
if [[ ! -s "$REGIONS" ]]; then
    echo "ERROR: region file not found or empty: $REGIONS" >&2
    exit 1
fi

find_bedtools
find_python

mkdir -p "$(dirname "$OUT_FASTA")"

TMP_BED=""
if [[ -n "$KEEP_BED" ]]; then
    mkdir -p "$(dirname "$KEEP_BED")"
    TMP_BED=$KEEP_BED
else
    TMP_BED=$(mktemp)
    trap 'rm -f "$TMP_BED"' EXIT
fi

"$PYTHON" - "$REGIONS" > "$TMP_BED" <<'PY'
import re
import sys
from collections import defaultdict

region_path = sys.argv[1]
region_re = re.compile(r"^([^:\s]+):([0-9]+)-([0-9]+)$")
regions = defaultdict(list)

with open(region_path) as handle:
    for line_no, line in enumerate(handle, 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        token = line.split()[0]
        match = region_re.match(token)
        if not match:
            raise SystemExit(
                "ERROR: {}:{} bad region {!r}; expected chr:start-end".format(
                    region_path, line_no, token
                )
            )
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        if start < 1 or end < 1:
            raise SystemExit("ERROR: {}:{} coordinates must be positive".format(region_path, line_no))
        if start > end:
            start, end = end, start
        regions[chrom].append((start - 1, end))

for chrom in sorted(regions):
    merged = []
    for start, end in sorted(regions[chrom]):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    for start, end in merged:
        print("{}\t{}\t{}".format(chrom, start, end))
PY

if [[ ! -s "$TMP_BED" ]]; then
    echo "ERROR: no valid regions were written to BED: $TMP_BED" >&2
    exit 1
fi

"$BEDTOOLS" maskfasta \
    -fi "$GENOME" \
    -bed "$TMP_BED" \
    -fo "$OUT_FASTA" \
    -mc "$MASK_CHAR"

echo "Genome:  $GENOME" >&2
echo "Regions: $REGIONS" >&2
echo "BED:     $TMP_BED" >&2
echo "Output:  $OUT_FASTA" >&2
