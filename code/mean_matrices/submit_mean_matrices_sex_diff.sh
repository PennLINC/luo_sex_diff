#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "HBN" "NKI")

atlases=("schaefer217")

r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/mean_matrices/mean_matrices_sex_diff.R"

# loop through each dataset and atlas
for dataset in "${datasets[@]}"; do    
    # where to save output and error logs
    logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/mean_matrices/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi
  
  for atlas in "${atlases[@]}"; do
      job_name="${dataset}_mean_matrices_sex_diff_${atlas}"
      
      # submit the job to SLURM using my Singularity container
      sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=24:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${atlas}"
      echo "Submitted ${dataset} ${atlas}"
  done
done
 