#!/bin/bash

# set variables
metrics=("GBC" "WNC" "BNC" "edge" "networkpair")
atlases=("schaefer200x17")

r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/fit_GAMs/fit_GAMs_all_datasets_clinical_maineffects_sexstratified.R"

for metric in "${metrics[@]}"; do
# where to save output and error logs
    logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/fit_GAMs/all_datasets/${metric}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi
    
    # loop through each metric
    for atlas in "${atlases[@]}"; do
        job_name="all_datasets_fit_GAMs_clinical_main_effects_${metric}_${atlas}"
        
        # submit the job to SLURM using my Singularity container
        sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=24:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${metric} ${atlas}"
        echo "Submitted ${dataset} ${metric}"
    done
done
 
 