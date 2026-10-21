# setting up the Conda Environment
# to create the `sex_diff_py311` environment, run:

module load miniconda3/25.5.1
cd ~/parcc/network_replication/covariate_analyses/sex_diff/software/python_env

conda env create \
    --prefix ~/parcc/network_replication/.conda/envs/sex_diff_py311 \
    -f sex_diff_py311_pinned.yml

# Tell Conda where the project's shared environments are located.
# This modifies the user's ~/.condarc and only needs to be done once.
conda config --add envs_dirs \
    ~/parcc/network_replication/.conda/envs

conda activate sex_diff_py311