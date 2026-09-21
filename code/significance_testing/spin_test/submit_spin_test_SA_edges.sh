#!/bin/bash

# set variables
datasets=("PNC" "HCPD" "NKI" "HBN" "all_datasets")
atlases=("schaefer200x17")
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/significance_testing/spin_test/spin_test_SA_edges.R"

# loop through each multisite dataset
for dataset in "${datasets[@]}"; do   
    for atlas in "${atlases[@]}"; do
        
        # where to save output and error logs
        logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/spin_test/${dataset}"
        if [ ! -d "${logs_dir}" ]; then
            mkdir -p ${logs_dir}
        fi

        # make outputs_root
        outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/${dataset}/spin_test"
        if [ ! -d "${outputs_root}" ]; then
            mkdir -p ${outputs_root}
        fi

        job_name="${dataset}_${atlas}_spin_test"
        
        # submit the job to SLURM using my Singularity container
        sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=12:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script} ${dataset} ${atlas}"
        echo "Submitted ${dataset} ${atlas} ${metric}"
   
  done
done