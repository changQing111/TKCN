#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
implementation="$project_root/estimate_copy_from_realcount_depth.sh"

if [[ ! -x "$implementation" ]]; then
    echo "ERROR: project implementation not found or not executable: $implementation" >&2
    exit 1
fi

exec "$implementation" "$@"
