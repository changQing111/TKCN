#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  build_copy_background.sh [-g genes.txt] [-d db_root] [-b base_db] [-o output_dir] [-t threads]

Build the fifth-step background for every target gene:
  target.non_target_copies.k31 = union of all copy KMC databases from other genes
  target.background.k31      = base background union target.non_target_copies.k31

Options:
  -g FILE   Gene list, one gene ID per line (default: multi_copy_gene.txt).
  -d DIR    Root containing <gene>/kmc_db/*.k31.kmc_pre/.kmc_suf
            (default: kmer_intersect).
  -b DB     Base-background KMC database prefix (default: union_regions.masked.k31).
  -o DIR    Output directory (default: copy_background_kmc).
  -t INT    Number of kmc_tools threads (default: 32).
  -h        Show this help message.

Environment:
  KMC_TOOLS_BIN  kmc_tools executable or path (default: kmc_tools).

Example:
  ./build_copy_background.sh \
      -g multi_copy_gene.txt \
      -d kmer_intersect \
      -b union_regions.masked.k31 \
      -o copy_background_kmc \
      -t 32
EOF
}

genes_file="multi_copy_gene.txt"
db_root="kmer_intersect"
base_db="union_regions.masked.k31"
output_dir="copy_background_kmc"
threads=32

while getopts ":g:d:b:o:t:h" opt; do
    case "$opt" in
        g) genes_file="$OPTARG" ;;
        d) db_root="$OPTARG" ;;
        b) base_db="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
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
            echo "Error: unknown option: -$OPTARG" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! -f "$genes_file" ]]; then
    echo "Error: gene list does not exist: $genes_file" >&2
    exit 1
fi

if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer: $threads" >&2
    exit 2
fi

kmc_tools_bin="${KMC_TOOLS_BIN:-kmc_tools}"
if ! command -v "$kmc_tools_bin" >/dev/null 2>&1; then
    echo "Error: kmc_tools was not found: $kmc_tools_bin" >&2
    echo "Set KMC_TOOLS_BIN to the kmc_tools executable path." >&2
    exit 127
fi

db_pair_exists() {
    local prefix="$1"
    [[ -f "${prefix}.kmc_pre" && -f "${prefix}.kmc_suf" ]]
}

require_db_pair() {
    local prefix="$1"
    if [[ ! -f "${prefix}.kmc_pre" || ! -f "${prefix}.kmc_suf" ]]; then
        echo "Error: incomplete KMC database pair: $prefix" >&2
        echo "Expected both ${prefix}.kmc_pre and ${prefix}.kmc_suf." >&2
        exit 1
    fi
}

copy_db_pair() {
    local source_prefix="$1"
    local target_prefix="$2"
    cp -- "${source_prefix}.kmc_pre" "${target_prefix}.kmc_pre"
    cp -- "${source_prefix}.kmc_suf" "${target_prefix}.kmc_suf"
}

move_db_pair() {
    local source_prefix="$1"
    local target_prefix="$2"
    mv -f -- "${source_prefix}.kmc_pre" "${target_prefix}.kmc_pre"
    mv -f -- "${source_prefix}.kmc_suf" "${target_prefix}.kmc_suf"
}

union_into() {
    local left_prefix="$1"
    local right_prefix="$2"
    local output_prefix="$3"

    "$kmc_tools_bin" -t"$threads" simple \
        "$left_prefix" \
        "$right_prefix" \
        union \
        "$output_prefix"
}

mapfile -t genes < <(awk 'NF && $1 !~ /^#/ {print $1}' "$genes_file")
if (( ${#genes[@]} < 2 )); then
    echo "Error: at least two genes are required in $genes_file." >&2
    exit 1
fi

declare -A seen_genes=()
for gene in "${genes[@]}"; do
    if [[ "$gene" == */* ]]; then
        echo "Error: gene IDs must not contain '/': $gene" >&2
        exit 1
    fi
    if [[ -n "${seen_genes[$gene]:-}" ]]; then
        echo "Error: duplicate gene ID in $genes_file: $gene" >&2
        exit 1
    fi
    seen_genes["$gene"]=1
done

require_db_pair "$base_db"
mkdir -p "$output_dir"

for target_gene in "${genes[@]}"; do
    non_target_prefix="$output_dir/${target_gene}.non_target_copies.k31"
    background_prefix="$output_dir/${target_gene}.background.k31"
    manifest="$output_dir/${target_gene}.non_target_copies.list"
    tmp_dir="$(mktemp -d "$output_dir/.${target_gene}.tmp.XXXXXX")"
    cleanup_tmp() {
        rm -rf -- "$tmp_dir"
    }
    trap cleanup_tmp EXIT

    for other_gene in "${genes[@]}"; do
        [[ "$other_gene" == "$target_gene" ]] && continue
        source_dir="$db_root/$other_gene/kmc_db"
        if [[ ! -d "$source_dir" ]]; then
            echo "Error: KMC directory does not exist: $source_dir" >&2
            exit 1
        fi
    done

    mapfile -t source_dbs < <(
        for other_gene in "${genes[@]}"; do
            [[ "$other_gene" == "$target_gene" ]] && continue
            source_dir="$db_root/$other_gene/kmc_db"
            find "$source_dir" -maxdepth 1 -type f -name '*.k31.kmc_pre' -print
        done | sort
    )

    if (( ${#source_dbs[@]} == 0 )); then
        echo "Error: no non-target copy KMC databases found for $target_gene." >&2
        exit 1
    fi

    source_prefixes=()
    for source_pre in "${source_dbs[@]}"; do
        source_prefix="${source_pre%.kmc_pre}"
        require_db_pair "$source_prefix"
        source_prefixes+=("$source_prefix")
    done

    printf '%s\n' "${source_prefixes[@]}" > "$manifest"

    work_prefix="$tmp_dir/non_target_copies.k31"
    copy_db_pair "${source_prefixes[0]}" "$work_prefix"
    for ((i = 1; i < ${#source_prefixes[@]}; i++)); do
        next_prefix="$tmp_dir/non_target_${i}.k31"
        union_into "$work_prefix" "${source_prefixes[$i]}" "$next_prefix"
        rm -f -- "${work_prefix}.kmc_pre" "${work_prefix}.kmc_suf"
        work_prefix="$next_prefix"
    done
    move_db_pair "$work_prefix" "$non_target_prefix"

    background_work_prefix="$tmp_dir/background.k31"
    union_into "$base_db" "$non_target_prefix" "$background_work_prefix"
    move_db_pair "$background_work_prefix" "$background_prefix"

    rm -rf -- "$tmp_dir"
    trap - EXIT
    echo "[$target_gene] non-target copies: ${#source_prefixes[@]} databases"
    echo "[$target_gene] background: $background_prefix"
done
