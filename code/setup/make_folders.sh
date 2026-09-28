#!/usr/bin/env bash

# Create directory structure needed for sex-differences replication.

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff"

mkdir -p \
    "$ROOT/input/PNC/sample_info" \
    "$ROOT/input/HCPD/sample_info" \
    "$ROOT/input/NKI/sample_info" \
    "$ROOT/input/HBN/sample_info" \
    "$ROOT/input/all_datasets/sample_info" \
    "$ROOT/output"

echo "Folder setup complete."