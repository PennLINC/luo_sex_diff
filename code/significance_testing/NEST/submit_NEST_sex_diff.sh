#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "NKI" "HBN")
datasets=("HBN")
atlases=("schaefer200x17")
metrics=("WNC")
networks_of_interest=("DefaultA" "LimbicA_TempPole")
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/significance_testing/NEST/NEST_sex_diff.R"

# loop through each multisite dataset
for dataset in "${datasets[@]}"; do   
    for atlas in "${atlases[@]}"; do
        for metric in "${metrics[@]}"; do
            for network_of_interest in "${networks_of_interest[@]}"; do
                # where to save output and error logs
                logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/NEST/${dataset}/${metric}"
                if [ ! -d "${logs_dir}" ]; then
                    mkdir -p ${logs_dir}
                fi

                # make outputs_root
                outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/${dataset}/NEST/${metric}"
                if [ ! -d "${outputs_root}" ]; then
                    mkdir -p ${outputs_root}
                fi

                job_name="${dataset}_NEST_${metric}_${atlas}_${network_of_interest}"
                
                # submit the job to SLURM using my Singularity container
                sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=48:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${metric} ${atlas} ${network_of_interest}"
                echo "Submitted ${dataset} ${atlas} ${metric}"
        done
    done
  done
done