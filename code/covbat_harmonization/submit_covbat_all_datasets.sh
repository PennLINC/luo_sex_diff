#!/usr/bin/env bash

# Submit CovBat harmonization for functional connectivity metrics
# from the pooled dataset (PNC, NKI, HCPD, HBN combined).

ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/covbat_harmonization/covbat_all_datasets.R"

atlases=("schaefer200x17")
metrics=("FC_strength" "BNC" "WNC" "edge" "networkpair")

for atlas in "${atlases[@]}"; do
    for metric in "${metrics[@]}"; do

        # Where to save output and error logs
        logs_dir="$PROJDIR/code/logs/covbat_harmonization/all_datasets/${metric}"
        mkdir -p "$logs_dir"

        # Make output directory
        outputs_root="$PROJDIR/output/all_datasets/${metric}"
        mkdir -p "$outputs_root"

        job_name="all_datasets_covbat_${metric}_${atlas}"

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
            --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$metric\" \"$atlas\""
        then
            echo "Submitted all_datasets ${atlas} ${metric} (${partition}, ${cpus} CPU(s), ${time_limit})"
        else
            echo "FAILED to submit all_datasets ${atlas} ${metric}" >&2
        fi

    done
done