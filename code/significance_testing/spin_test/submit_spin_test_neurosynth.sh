#!/bin/bash

# set variables
r_script="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/significance_testing/spin_test/spin_test_neurosynth.R"


# where to save output and error logs
logs_dir="/cbica/projects/network_replication/covariate_analyses/sex_diff/code/logs/spin_test/all_datasets"
if [ ! -d "${logs_dir}" ]; then
    mkdir -p ${logs_dir}
fi

# make outputs_root
outputs_root="/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets/spin_test"
if [ ! -d "${outputs_root}" ]; then
    mkdir -p ${outputs_root}
fi

job_name="neurosynth_spin_test"

# submit the job to SLURM using my Singularity container
sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=12:00:00 --propagate=NONE --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/network_replication/covariate_analyses/sex_diff/software/r_packages/r-packages-for-cubic_0.1.0.sif Rscript --save ${r_script}"
echo "Submitted neurosynth spin test"

