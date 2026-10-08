


#pull docker image with necessary R packages
#docker site: https://hub.docker.com/r/audreycluo/r-packages-for-parcc
#tag: https://hub.docker.com/layers/audreycluo/r-packages-for-parcc/0.0.1/images/sha256-b59c4703bfce53aaf6d685d139c966d916cda5a573b5b021368434e45cd43d80

#pull the docker image onto parcc cluster:
module load apptainer/1.4.1 # if on parcc
cd /ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/software
mkdir -p r_packages
singularity pull r_packages/r-packages-for-parcc_0.0.1.sif docker://audreycluo/r-packages-for-parcc:0.0.1
