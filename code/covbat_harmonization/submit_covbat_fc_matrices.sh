#!/bin/bash

# set variables
datasets=("HCPD" "HBN" "all_datasets")
atlases=("schaefer200x17")
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/covbat_harmonization/covbat_fc_matrices.R"

mem="32G"
time_limit="24:00:00"

# loop through each dataset (single-site datasets don't need this script --
# only include ones with a site variable: HCPD, HBN, or the pooled run)
for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do
        # where to save output and error logs
        logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/covbat_harmonization/${dataset}/fc_matrices"
        if [ ! -d "${logs_dir}" ]; then
            mkdir -p "${logs_dir}"
        fi

        # make outputs_root
        outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/${dataset}/fc_matrices"
        if [ ! -d "${outputs_root}" ]; then
            mkdir -p "${outputs_root}"
        fi

        job_name="${dataset}_covbat_fc_matrices_${atlas}"

        # submit the job to SLURM using the Singularity container
        sbatch --job-name=${job_name} \
            --nodes=1 --ntasks=1 --cpus-per-task=4 \
            --mem=${mem} \
            --time=${time_limit} \
            --propagate=NONE \
            --output=${logs_dir}/${job_name}_%j.out \
            --error=${logs_dir}/${job_name}_%j.err \
            --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${atlas}"

        echo "Submitted ${dataset} ${atlas} fc_matrices with ${mem} memory"
    done
done
