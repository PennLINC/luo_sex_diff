#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "HBN" "NKI")
atlases=("schaefer200x7" "schaefer200x17")
metrics=("GBC" "BNC" "WNC" "edge" "networkpair")

r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/connectivity_measures/compute_connectivity_metrics.R"

# loop through each dataset and num_networks
for dataset in "${datasets[@]}"; do   
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do
            # where to save output and error logs
            logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/connectivity_measures/${dataset}/${metric}"
            if [ ! -d "${logs_dir}" ]; then
                mkdir -p ${logs_dir}
            fi

            # make outputs_root
            outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/${dataset}/${metric}"
            if [ ! -d "${outputs_root}" ]; then
                mkdir -p ${outputs_root}
            fi

            job_name="${dataset}_compute${metric}_${atlas}"
            
            # submit the job to SLURM using my Singularity container
            sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=24:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${metric} ${atlas}"
            echo "Submitted ${dataset} ${atlas} ${metric}"
    done
  done
done
 