#!/usr/bin/env bash

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"
APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/connectivity_measures/compute_connectivity_metrics.R"

datasets=("PNC" "HCPD" "HBN" "NKI")
atlases=("schaefer200x17")
metrics=("FC_strength" "BNC" "WNC" "edge" "networkpair")


for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do

            logs_dir="$PROJDIR/code/logs/connectivity_measures/${dataset}/${metric}"
            mkdir -p "$logs_dir"

            outputs_root="$PROJDIR/output/${dataset}/${metric}"
            mkdir -p "$outputs_root"

            job_name="${dataset}_compute${metric}_${atlas}"

            if sbatch \
                --job-name="$job_name" \
                --partition=genoa-std-mem \
                --nodes=1 \
                --ntasks=1 \
                --cpus-per-task=1 \
                --time=01:00:00 \
                --output="${logs_dir}/${job_name}_%j.out" \
                --error="${logs_dir}/${job_name}_%j.err" \
                --wrap="$APPTAINER exec --cleanenv --bind ${ROOT}:${ROOT} ${SIF} Rscript ${r_script} ${dataset} ${metric} ${atlas}"
            then
                echo "Submitted ${dataset} ${atlas} ${metric}"
            else
                echo "FAILED to submit ${dataset} ${atlas} ${metric}" >&2
            fi
        done
    done
done