#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "HBN" "NKI" "all_datasets")
metrics=("GBC" "WNC" "BNC" "edge" "networkpair")
atlases=("schaefer200x17")

r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/compute_effect_size/resi_sex_diff.R"

# loop through each dataset and metric
for dataset in "${datasets[@]}"; do    
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do
            # where to save output and error logs
            logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/compute_effect_size/${dataset}/${metric}"
            if [ ! -d "${logs_dir}" ]; then
                mkdir -p ${logs_dir}
            fi
        
            
            job_name="${dataset}_compute_effect_size_${metric}_${atlas}"
            
          
            # submit the job to SLURM
            sbatch --job-name=${job_name} \
                   --nodes=1 --ntasks=1 --cpus-per-task=2 \
                   --mem=3G --time=2:00:00 --propagate=NONE \
                   --output=${logs_dir}/${job_name}_%j.out \
                   --error=${logs_dir}/${job_name}_%j.err \
                   --wrap="singularity run --cleanenv \
                          /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.1.sif \
                          Rscript --save ${r_script} ${dataset} ${metric} ${atlas}"
            echo "Submitted ${dataset} ${metric} ${atlas}"
    done
  done
done
 