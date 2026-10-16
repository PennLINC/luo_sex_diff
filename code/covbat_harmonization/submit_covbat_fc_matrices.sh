#!/usr/bin/env bash

# submit covbat for connectivity matrices in HBN, HCP-D, and pooled dataset

ROOT="~/parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/covbat_harmonization/covbat_fc_matrices.R"

datasets=("HCPD" "HBN" "all_datasets")
atlases=("schaefer200x17")

partition="genoa-lrg-mem"
cpus=2
time_limit="02:00:00"

# Loop through each dataset.
# Single-site datasets don't need this script; only datasets with a site
# variable are included (HCPD, HBN, and the pooled all_datasets run).
for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do

        # Where to save output and error logs
        logs_dir="$PROJDIR/code/logs/covbat_harmonization/${dataset}/fc_matrices"
        mkdir -p "$logs_dir"

        # Make output directory
        outputs_root="$PROJDIR/output/${dataset}/fc_matrices"
        mkdir -p "$outputs_root"

        job_name="${dataset}_covbat_fc_matrices_${atlas}"

        if sbatch \
            --job-name="$job_name" \
            --partition="$partition" \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="$cpus" \
            --time="$time_limit" \
            --output="${logs_dir}/${job_name}_%j.out" \
            --error="${logs_dir}/${job_name}_%j.err" \
            --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$dataset\" \"$atlas\""
        then
            echo "Submitted ${dataset} ${atlas} fc_matrices (${partition}, ${cpus} CPUs)"
        else
            echo "FAILED to submit ${dataset} ${atlas} fc_matrices" >&2
        fi

    done
done