#!/usr/bin/env bash
# run AFTER group_template.py has completed for each dataset

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

PYTHON="$ROOT/.conda/envs/sex_diff_py311/bin/python"
py_script="$PROJDIR/code/gradient_analysis/procrustes_alignment.py"

datasets=("PNC" "NKI" "HCPD" "HBN" "all_datasets")

cpus=1
time_limit="00:10:00"

for dataset in "${datasets[@]}"; do
    logs_dir="$PROJDIR/code/logs/gradient_analysis/${dataset}"
    mkdir -p "$logs_dir"

    job_name="procrustes_alignment_${dataset}"

    if sbatch \
        --job-name="$job_name" \
        --partition=genoa-std-mem \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$cpus" \
        --time="$time_limit" \
        --output="${logs_dir}/${job_name}_%j.out" \
        --error="${logs_dir}/${job_name}_%j.err" \
        --wrap="\"$PYTHON\" -u \"$py_script\" \"$dataset\""
    then
        echo "Submitted ${job_name} (${cpus} CPU(s), ${time_limit})"
    else
        echo "FAILED to submit ${job_name}" >&2
    fi
done