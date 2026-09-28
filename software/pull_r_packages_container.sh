


#pull docker image with necessary R packages
#docker site: https://hub.docker.com/r/audreycluo/r-packages-for-cubic/tags 
#tag: https://hub.docker.com/layers/audreycluo/r-packages-for-cubic/0.1.1/images/sha256-2209677ec5857243d6411e87171c6e9c1be02209a756efe3ca669ea9be70a1b8

#pull the docker image onto cubic cluster:
cd /cbica/projects/network_replication/covariate_analyses/sex_diff/software
mkdir -p r_packages
singularity pull r_packages/r-packages-for-cubic_0.1.1.sif docker://audreycluo/r-packages-for-cubic:0.1.1
