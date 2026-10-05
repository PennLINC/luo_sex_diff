#!/usr/bin/env bash

# Submit Schaefer200x17 parcellation rotations for spin testing
ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/spin_test/rotate_schaefer200x17.R"

# Where to save output and error logs
logs_dir="$PROJDIR/code/logs/spin_test/rotate_parcellation"
mkdir -p "$logs_dir"

job_name="rotate_parcellation_schaefer200x17"
cpus=1
time_limit="00:30:00"

if sbatch \
    --job-name="$job_name" \
    --partition=genoa-std-mem \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$cpus" \
    --time="$time_limit" \
    --output="${logs_dir}/${job_name}_%j.out" \
    --error="${logs_dir}/${job_name}_%j.err" \
    --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\""
then
    echo "Submitted ${job_name} (${cpus} CPU(s), ${time_limit})"
else
    echo "FAILED to submit ${job_name}" >&2
fi
 