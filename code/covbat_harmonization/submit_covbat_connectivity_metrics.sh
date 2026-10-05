#!/usr/bin/env bash

# Submit CovBat for functional connectivity metrics in HBN and HCP-D,
# the multi-site datasets.

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/covbat_harmonization/covbat_connectivity_metrics.R"

datasets=("HCPD" "HBN")
atlases=("schaefer200x17")
metrics=("FC_strength" "BNC" "WNC" "edge" "networkpair")

# Loop through each multisite dataset, atlas, and metric
for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do

            # Where to save output and error logs
            logs_dir="$PROJDIR/code/logs/covbat_harmonization/${dataset}/${metric}"
            mkdir -p "$logs_dir"

            # Make output directory
            outputs_root="$PROJDIR/output/${dataset}/${metric}"
            mkdir -p "$outputs_root"

            job_name="${dataset}_covbat_${metric}_${atlas}"

            # Edge requires more memory; otherwise minimize resources requested
            if [ "$metric" = "edge" ]; then
                partition="genoa-lrg-mem"
                cpus=2
                time_limit="01:00:00"
            else
                partition="genoa-std-mem"
                cpus=1
                time_limit="00:30:00"
            fi

            if sbatch \
                --job-name="$job_name" \
                --partition="$partition" \
                --nodes=1 \
                --ntasks=1 \
                --cpus-per-task="$cpus" \
                --time="$time_limit" \
                --output="${logs_dir}/${job_name}_%j.out" \
                --error="${logs_dir}/${job_name}_%j.err" \
                --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$dataset\" \"$metric\" \"$atlas\""
            then
                echo "Submitted ${dataset} ${atlas} ${metric} (${partition}, ${cpus} CPU(s), ${time_limit})"
            else
                echo "FAILED to submit ${dataset} ${atlas} ${metric}" >&2
            fi

        done
    done
done