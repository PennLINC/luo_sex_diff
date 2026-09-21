#!/bin/bash

# set variables
atlases=("schaefer200x7" "schaefer200x17")
metrics=("GBC" "BNC" "WNC" "edge" "networkpair")
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/covbat_harmonization/covbat_all_datasets.R"

# loop through each atlas/metric
for atlas in "${atlases[@]}"; do
    for metric in "${metrics[@]}"; do
    # where to save output and error logs
    logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/covbat_harmonization/all_datasets/${metric}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p "${logs_dir}"
    fi

    # make outputs_root
    outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets/${metric}"
    if [ ! -d "${outputs_root}" ]; then
        mkdir -p "${outputs_root}"
    fi

    job_name="all_datasets_covbat_${metric}_${atlas}"

    # set memory conditionally
    if [ "${metric}" == "edge" ]; then
        mem="32G"
    else
        mem="8G"
    fi

    # submit the job to SLURM
    sbatch --job-name=${job_name} \
        --nodes=1 --ntasks=1 --cpus-per-task=4 \
        --mem=${mem} \
        --time=24:00:00 \
        --propagate=NONE \
        --output=${logs_dir}/${job_name}_%j.out \
        --error=${logs_dir}/${job_name}_%j.err \
        --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${metric} ${atlas}"

    echo "Submitted all datasets ${atlas} ${metric} with ${mem} memory"
done

done
