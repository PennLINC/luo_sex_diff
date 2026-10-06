#!/usr/bin/env bash
# run AFTER compute_dispersion.py has completed for each dataset

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/gradient_analysis/resi_dispersion.R"

datasets=("PNC" "NKI" "HCPD" "HBN" "all_datasets")

cpus=2
time_limit="00:30:00"

for dataset in "${datasets[@]}"; do
    logs_dir="$PROJDIR/code/logs/gradient_analysis/${dataset}"
    mkdir -p "$logs_dir"

    outputs_root="$PROJDIR/output/${dataset}/gradient_dispersion/resi"
    mkdir -p "$outputs_root"

    job_name="resi_dispersion_${dataset}"

    if sbatch \
        --job-name="$job_name" \
        --partition=genoa-std-mem \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$cpus" \
        --time="$time_limit" \
        --output="${logs_dir}/${job_name}_%j.out" \
        --error="${logs_dir}/${job_name}_%j.err" \
        --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$dataset\""
    then
        echo "Submitted ${job_name} (${cpus} CPU(s), ${time_limit})"
    else
        echo "FAILED to submit ${job_name}" >&2
    fi
 done