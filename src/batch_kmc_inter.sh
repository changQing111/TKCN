#!/usr/bin/env bash
set -euo pipefail

kmc_dir="TaFT-B1_all_copy_single_kmc"
out="TaFT-B1_all_copy_single_kmc.ALL_INTERSECT.k31"
KMC_TOOLS=${KMC_TOOLS:-kmc_tools}

usage() {
  cat <<EOF
Usage:
  bash $0 [options] [INPUT_DIR OUTPUT_PREFIX]

Options:
  -i, --input-dir DIR   Directory containing *.kmc_pre/*.kmc_suf databases.
                        Default: $kmc_dir
  -o, --output PREFIX   Output KMC database prefix.
                        Default: $out
  --kmc-tools PATH      kmc_tools executable. Default: $KMC_TOOLS
  -h, --help            Show this help.

Examples:
  bash $0 kmc_db kmc_db/target_shared.k31
  bash $0 --input-dir kmc_db --output kmc_db/target_shared.k31
EOF
}

positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a directory" >&2; exit 1; }
      kmc_dir=$2
      shift 2
      ;;
    --input-dir=*)
      kmc_dir=${1#*=}
      shift
      ;;
    -o|--output)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an output prefix" >&2; exit 1; }
      out=$2
      shift 2
      ;;
    --output=*)
      out=${1#*=}
      shift
      ;;
    --kmc-tools)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an executable path" >&2; exit 1; }
      KMC_TOOLS=$2
      shift 2
      ;;
    --kmc-tools=*)
      KMC_TOOLS=${1#*=}
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

if [[ ${#positional[@]} -gt 2 ]]; then
  echo "ERROR: too many positional arguments" >&2
  usage >&2
  exit 1
fi

if [[ ${#positional[@]} -ge 1 ]]; then
  kmc_dir=${positional[0]}
fi

if [[ ${#positional[@]} -ge 2 ]]; then
  out=${positional[1]}
fi

if [[ ! -d "$kmc_dir" ]]; then
  echo "ERROR: input directory not found: $kmc_dir" >&2
  exit 1
fi

if ! command -v "$KMC_TOOLS" >/dev/null 2>&1; then
  echo "ERROR: kmc_tools executable not found: $KMC_TOOLS" >&2
  exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

mkdir -p "$(dirname "$out")"
out_abs=$(cd "$(dirname "$out")" && pwd -P)/$(basename "$out")

mapfile -t found_dbs < <(
  find "$kmc_dir" -maxdepth 1 -name "*.kmc_pre" \
    | sed 's/\.kmc_pre$//' \
    | sort
)

dbs=()
for db in "${found_dbs[@]}"; do
  db_abs=$(cd "$(dirname "$db")" && pwd -P)/$(basename "$db")
  if [[ "$db_abs" == "$out_abs" ]]; then
    continue
  fi
  dbs+=("$db")
done

if [[ ${#dbs[@]} -eq 0 ]]; then
  echo "ERROR: no KMC databases found in $kmc_dir" >&2
  exit 1
fi

for db in "${dbs[@]}"; do
  if [[ ! -s "${db}.kmc_suf" ]]; then
    echo "ERROR: missing KMC suffix file for ${db}" >&2
    exit 1
  fi
done

acc="${dbs[0]}"

for i in $(seq 1 $((${#dbs[@]} - 1))); do
  next="$tmpdir/intersect_${i}.k31"
  "$KMC_TOOLS" simple "$acc" "${dbs[$i]}" intersect "$next"
  acc="$next"
done

cp "${acc}.kmc_pre" "${out}.kmc_pre"
cp "${acc}.kmc_suf" "${out}.kmc_suf"

"$KMC_TOOLS" transform "$out" dump -s "${out}.tsv"

echo "Intersect KMC database written to: $out" >&2
echo "Intersect TSV written to: ${out}.tsv" >&2
