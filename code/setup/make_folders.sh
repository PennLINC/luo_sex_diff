#!/usr/bin/env bash
# Create directory structure needed for sex-differences project.

set -euo pipefail  # stop on the first error

ROOT="$HOME/parcc/network_replication/covariate_analyses/sex_diff"

dirs=(
    "$ROOT/input/PNC/sample_info"
    "$ROOT/input/HCPD/sample_info"
    "$ROOT/input/NKI/sample_info"
    "$ROOT/input/HBN/sample_info"
    "$ROOT/input/all_datasets/sample_info"
    "$ROOT/output"
)

mkdir -p "${dirs[@]}"

# confirm every folder exists before reporting success
for d in "${dirs[@]}"; do
    if [ ! -d "$d" ]; then
        echo "ERROR: failed to create $d" >&2
        exit 1
    fi
done

echo "Folder setup complete under $ROOT"