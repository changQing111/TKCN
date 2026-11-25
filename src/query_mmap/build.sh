#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
out_bin=${1:-"$script_dir/query_mmap"}

pick_cxx() {
  if [[ -n "${CXX:-}" ]]; then
    printf '%s' "$CXX"
    return
  fi
  for candidate in \
    /data/changq/miniforge3/bin/g++ \
    /data/changq/miniforge3/bin/x86_64-conda-linux-gnu-g++ \
    g++ \
    c++; do
    if command -v "$candidate" >/dev/null 2>&1; then
      if "$candidate" -x c++ -std=c++14 -fsyntax-only - <<< '#include <memory>
int main(){ auto p = std::make_unique<int>(0); return *p; }' 2>/dev/null; then
        printf '%s' "$candidate"
        return
      fi
    fi
  done
  echo "ERROR: C++14-capable compiler not found; set CXX=/path/to/g++" >&2
  exit 1
}

cxx=$(pick_cxx)

"$python_bin" "$script_dir/build_mmap_api.py"
ln -sfn mmap_api "$script_dir/kmc_api"

extra_flags=()
if "$cxx" -x c++ -std=c++14 -static-libstdc++ -static-libgcc -o /dev/null - <<< 'int main(){return 0;}' 2>/dev/null; then
  extra_flags=(-static-libstdc++ -static-libgcc)
fi

"$cxx" -O3 -std=c++14 \
  "${extra_flags[@]}" \
  -I "$script_dir/mmap_api" \
  "$script_dir/query_mmap.cpp" \
  "$script_dir/mmap_api/kmc_file.cpp" \
  "$script_dir/mmap_api/kmer_api.cpp" \
  "$script_dir/mmap_api/mmer.cpp" \
  -o "$out_bin"

echo "Built query_mmap: $out_bin (compiler: $cxx)" >&2
