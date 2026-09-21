#!/bin/bash
# submit with ./submit_compute_dispersion.sh
# run AFTER procrustes_alignment.py has completed for each dataset

datasets=("PNC" "NKI" "HCPD" "HBN" "all_datasets")
script_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/gradient_analysis"

source /cbica/projects/network_replication/miniconda3/etc/profile.d/conda.sh
conda activate sex_diff_py311

for dataset in "${datasets[@]}"; do
    logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/gradient_analysis/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p "${logs_dir}"
    fi

    job_name="compute_dispersion_${dataset}"

    sbatch --job-name=${job_name} \
        --nodes=1 --ntasks=1 --cpus-per-task=4 \
        --mem=16G \
        --time=4:00:00 \
        --output=${logs_dir}/${job_name}_%j.out \
        --error=${logs_dir}/${job_name}_%j.err \
        --wrap="python ${script_dir}/compute_dispersion.py ${dataset}"

    echo "Submitted compute_dispersion for ${dataset}"
done
