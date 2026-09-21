#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "HBN" "NKI")
metrics=("GBC" "WNC" "BNC" "edge" "networkpair")
atlases=("schaefer200x17")
#atlases=("schaefer200x7" "schaefer200x17")

r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/fit_GAMs/fit_GAMs_sex_diff.R"

# loop through each dataset and metric
for dataset in "${datasets[@]}"; do    
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do
            # where to save output and error logs
            logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/fit_GAMs/${dataset}/${metric}"
            if [ ! -d "${logs_dir}" ]; then
                mkdir -p ${logs_dir}
            fi
        
            
            job_name="${dataset}_fit_GAMs_sex_diff_${metric}_${atlas}"
            
            # submit the job to SLURM using my Singularity container
            sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=24:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${metric} ${atlas}"
            echo "Submitted ${dataset} ${metric} ${atlas}"
    done
  done
done
 