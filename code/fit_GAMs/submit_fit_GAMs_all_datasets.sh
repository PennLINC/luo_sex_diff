#!/usr/bin/env bash

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/fit_GAMs/fit_GAMs_all_datasets.R"
metrics=("FC_strength" "WNC" "BNC" "edge" "networkpair")
atlases=("schaefer200x17")

for metric in "${metrics[@]}"; do

    logs_dir="$PROJDIR/code/logs/fit_GAMs/all_datasets/${metric}"
    mkdir -p "$logs_dir"

    for atlas in "${atlases[@]}"; do

        job_name="all_datasets_fit_GAMs_${metric}_${atlas}"


        # Edge needs more resources
        if [ "$metric" = "edge" ]; then
            cpus=4
            time_limit="02:00:00"
        else
            cpus=4
            time_limit="00:30:00"   
        fi
        
        if sbatch \
            --job-name="$job_name" \
            --partition=genoa-std-mem \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="$cpus" \
            --time="$time_limit" \
            --output="${logs_dir}/${job_name}_%j.out" \
            --error="${logs_dir}/${job_name}_%j.err" \
            --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$metric\" \"$atlas\""
        then
            echo "Submitted all_datasets ${metric} ${atlas} (${cpus} CPUs, ${time_limit})"
        else
            echo "FAILED to submit all_datasets ${metric} ${atlas}" >&2
        fi

    done
done