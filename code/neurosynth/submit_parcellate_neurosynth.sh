#!/usr/bin/env bash

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

PYTHON="$ROOT/.conda/envs/sex_diff_py311/bin/python"
py_script="$PROJDIR/code/neurosynth/parcellate_neurosynth.py"

logs_dir="$PROJDIR/code/logs/neurosynth"
mkdir -p "$logs_dir"

job_name="parcellate_neurosynth"

cpus=2
time_limit="08:00:00"

if sbatch \
    --job-name="$job_name" \
    --partition=genoa-std-mem \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$cpus" \
    --time="$time_limit" \
    --output="${logs_dir}/${job_name}_%j.out" \
    --error="${logs_dir}/${job_name}_%j.err" \
    --wrap="\"$PYTHON\" -u \"$py_script\""
then
    echo "Submitted ${job_name} (${cpus} CPU(s), ${time_limit})"
else
    echo "FAILED to submit ${job_name}" >&2
fi