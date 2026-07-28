#!/usr/bin/env bash
set -euo pipefail

MINIPROT_BIN=${MINIPROT_BIN:-miniprot}
GENOME_DIR=${GENOME_DIR:-/data/genome}
THREADS=${THREADS:-16}
QUERY=""
VARIETY=""
GENOME=""
INDEX=""
GFF=""

usage() {
    cat <<'EOF'
Usage:
  run_miniport.sh -q QUERY.pep.fa -v VARIETY [options]

Options:
  -q FILE   Protein query FASTA.
  -v NAME   Variety name. Used with GENOME_DIR/NAME.fa by default.
  -g FILE   Genome FASTA. Overrides GENOME_DIR/NAME.fa.
  -i FILE   Miniprot index prefix. Default: VARIETY_genome.mpi.
  -o FILE   GFF output path. Default: TaFT-B1_VARIETY.gff.
  -t INT    Number of threads. Default: 16.
  --miniprot PATH
            Miniprot executable. Default: MINIPROT_BIN or miniprot in PATH.
  -h        Show this help message.

Environment:
  GENOME_DIR=/path/to/genomes
  MINIPROT_BIN=/path/to/miniprot

Example:
  run_miniport.sh \
    -q test/TaFT-B1_pep.fasta \
    -v CS \
    -g /data/genome/CS.fa \
    -o result/TaFT-B1_CS.gff
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q|--query)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a file" >&2; exit 2; }
            QUERY=$2
            shift 2
            ;;
        -v|--variety)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a name" >&2; exit 2; }
            VARIETY=$2
            shift 2
            ;;
        -g|--genome)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a FASTA" >&2; exit 2; }
            GENOME=$2
            shift 2
            ;;
        -i|--index)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a prefix" >&2; exit 2; }
            INDEX=$2
            shift 2
            ;;
        -o|--gff)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an output path" >&2; exit 2; }
            GFF=$2
            shift 2
            ;;
        -t|--threads)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an integer" >&2; exit 2; }
            THREADS=$2
            shift 2
            ;;
        --miniprot)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an executable" >&2; exit 2; }
            MINIPROT_BIN=$2
            shift 2
            ;;
        --miniprot=*)
            MINIPROT_BIN=${1#*=}
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$QUERY" || -z "$VARIETY" ]]; then
    echo "ERROR: -q/--query and -v/--variety are required" >&2
    usage >&2
    exit 2
fi

GENOME=${GENOME:-"${GENOME_DIR}/${VARIETY}.fa"}
INDEX=${INDEX:-"${VARIETY}_genome.mpi"}
GFF=${GFF:-"TaFT-B1_${VARIETY}.gff"}

[[ -f "$QUERY" ]] || { echo "ERROR: query FASTA not found: $QUERY" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: genome FASTA not found: $GENOME" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be a positive integer" >&2; exit 2; }

if [[ "$MINIPROT_BIN" == */* ]]; then
    [[ -x "$MINIPROT_BIN" ]] || { echo "ERROR: miniprot is not executable: $MINIPROT_BIN" >&2; exit 1; }
else
    command -v "$MINIPROT_BIN" >/dev/null 2>&1 || {
        echo "ERROR: miniprot executable not found: $MINIPROT_BIN" >&2
        echo "Set MINIPROT_BIN=/path/to/miniprot or use --miniprot." >&2
        exit 127
    }
fi

mkdir -p "$(dirname "$INDEX")" "$(dirname "$GFF")"

"$MINIPROT_BIN" -t "$THREADS" -d "$INDEX" "$GENOME"
"$MINIPROT_BIN" -I -u -t "$THREADS" --gff "$INDEX" "$QUERY" > "$GFF"

echo "Miniprot index written with prefix: $INDEX" >&2
echo "Miniprot GFF written to: $GFF" >&2
