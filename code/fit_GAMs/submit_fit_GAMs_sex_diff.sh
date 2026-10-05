#!/usr/bin/env bash

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/fit_GAMs/fit_GAMs_sex_diff.R"

datasets=("PNC" "HCPD" "HBN" "NKI")
metrics=("FC_strength" "WNC" "BNC" "edge" "networkpair")
atlases=("schaefer200x17")

# Loop through each dataset, atlas, and metric
for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do

            # Where to save output and error logs
            logs_dir="$PROJDIR/code/logs/fit_GAMs/${dataset}/${metric}"
            mkdir -p "$logs_dir"

            job_name="${dataset}_fit_GAMs_sex_diff_${metric}_${atlas}"

            # Edge GAMs take substantially longer than the other metrics.
            if [ "$metric" = "edge" ]; then
                time_limit="01:30:00"
            else
                time_limit="00:30:00"
            fi

            # Submit job to Slurm using the Apptainer R container
            if sbatch \
                --job-name="$job_name" \
                --partition=genoa-std-mem \
                --nodes=1 \
                --ntasks=1 \
                --cpus-per-task=4 \
                --time="$time_limit" \
                --output="${logs_dir}/${job_name}_%j.out" \
                --error="${logs_dir}/${job_name}_%j.err" \
                --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$dataset\" \"$metric\" \"$atlas\""
            then
                echo "Submitted ${dataset} ${metric} ${atlas} (4 CPUs, ${time_limit})"
            else
                echo "FAILED to submit ${dataset} ${metric} ${atlas}" >&2
            fi

        done
    done
done