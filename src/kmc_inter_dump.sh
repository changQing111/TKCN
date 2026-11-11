#!/usr/bin/env bash
set -euo pipefail

KMC_TOOLS=${KMC_TOOLS:-${KMC_TOOLS_BIN:-/data/changq/bio_soft/KMC3.2.4/bin/kmc_tools}}
out_tsv=""
force=0

usage() {
  cat <<EOF
Usage:
  bash $0 [options] SET1 SET2 OUTPUT_PREFIX

Arguments:
  SET1             First KMC database prefix to intersect.
  SET2             Second KMC database prefix to intersect.
  OUTPUT_PREFIX    Output KMC database prefix.

Options:
  --tsv PATH        Output TSV path from kmc_tools transform dump.
                    Default: OUTPUT_PREFIX.tsv
  --kmc-tools PATH  kmc_tools executable. Default: $KMC_TOOLS
  -f, --force       Overwrite existing output files.
  -h, --help        Show this help.

Examples:
  bash $0 sample.k31 masked_union.k31 sample_specific.k31
  bash $0 --tsv tables/sample_specific.tsv sample.k31 masked_union.k31 out/sample_specific.k31
EOF
}

positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tsv)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a path" >&2; exit 1; }
      out_tsv=$2
      shift 2
      ;;
    --tsv=*)
      out_tsv=${1#*=}
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
    -f|--force)
      force=1
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

if [[ ${#positional[@]} -ne 3 ]]; then
  echo "ERROR: SET1, SET2, and OUTPUT_PREFIX are required" >&2
  usage >&2
  exit 1
fi

set1=${positional[0]}
set2=${positional[1]}
out=${positional[2]}

if [[ -z "$out_tsv" ]]; then
  out_tsv="${out}.tsv"
fi

if ! command -v "$KMC_TOOLS" >/dev/null 2>&1; then
  echo "ERROR: kmc_tools executable not found: $KMC_TOOLS" >&2
  exit 1
fi

for db in "$set1" "$set2"; do
  if [[ ! -s "${db}.kmc_pre" ]]; then
    echo "ERROR: missing KMC prefix file: ${db}.kmc_pre" >&2
    exit 1
  fi
  if [[ ! -s "${db}.kmc_suf" ]]; then
    echo "ERROR: missing KMC suffix file: ${db}.kmc_suf" >&2
    exit 1
  fi
done

output_paths=(
  "$(readlink -m "${out}.kmc_pre")"
  "$(readlink -m "${out}.kmc_suf")"
  "$(readlink -m "$out_tsv")"
)
if [[ "${output_paths[2]}" == "${output_paths[0]}" || "${output_paths[2]}" == "${output_paths[1]}" ]]; then
  echo "ERROR: output TSV must differ from output database files: $out_tsv" >&2
  exit 1
fi
for db in "$set1" "$set2"; do
  input_paths=(
    "$(readlink -m "${db}.kmc_pre")"
    "$(readlink -m "${db}.kmc_suf")"
  )
  for output_path in "${output_paths[@]}"; do
    for input_path in "${input_paths[@]}"; do
      if [[ "$output_path" == "$input_path" ]]; then
        echo "ERROR: output file collides with input database file: $output_path" >&2
        exit 1
      fi
    done
  done
done

mkdir -p "$(dirname "$out")" "$(dirname "$out_tsv")"

if [[ $force -eq 0 ]] && {
  [[ -e "${out}.kmc_pre" ]] || [[ -e "${out}.kmc_suf" ]] || [[ -e "$out_tsv" ]]
}; then
  echo "ERROR: output already exists: ${out}.kmc_pre, ${out}.kmc_suf, or $out_tsv" >&2
  echo "Use --force to overwrite." >&2
  exit 1
fi

if [[ $force -eq 1 ]]; then
  rm -f "${out}.kmc_pre" "${out}.kmc_suf" "$out_tsv"
fi

"$KMC_TOOLS" simple \
  "$set1" \
  "$set2" \
  intersect \
  "$out"

"$KMC_TOOLS" transform \
  "$out" \
  dump \
  -s \
  "$out_tsv"

echo "Intersected KMC database written to: $out" >&2
echo "Intersected TSV written to: $out_tsv" >&2
