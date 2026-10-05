#!/bin/bash
# submit with ./submit_resi_dispersion.sh
# run AFTER compute_dispersion.py has completed for each dataset

datasets=("PNC" "NKI" "HCPD" "HBN" "all_datasets")
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/gradient_analysis/resi_dispersion.R"

for dataset in "${datasets[@]}"; do
    logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/gradient_analysis/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p "${logs_dir}"
    fi

    outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/${dataset}/gradient_dispersion/resi"
    if [ ! -d "${outputs_root}" ]; then
        mkdir -p "${outputs_root}"
    fi

    job_name="resi_dispersion_${dataset}"

    sbatch --job-name=${job_name} \
        --nodes=1 --ntasks=1 --cpus-per-task=4 \
        --mem=8G \
        --time=4:00:00 \
        --output=${logs_dir}/${job_name}_%j.out \
        --error=${logs_dir}/${job_name}_%j.err \
        --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-parcc_0.0.1.sif Rscript --save ${r_script} ${dataset}"

    echo "Submitted resi_dispersion for ${dataset}"
done
