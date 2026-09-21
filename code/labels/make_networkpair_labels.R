################## 
# Set Variables 
################## 
dataset = "PNC" # any dataset works
measure = "networkpair"
 
################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
outputs_root <- config_data$covariate_output_root
measure_outputs_dir <- paste0(outputs_root, "/", measure)
 
################## 
# Write out labels
################## 
num_networks = 7
schaefer200x7_networkpair <- read.csv(sprintf("%s/networkpair_subxnetpair_matrix_schaefer200x%s_orig.csv", measure_outputs_dir, num_networks))
schaefer200x7_networkpair_labels <- setdiff(names(schaefer200x7_networkpair), "subject")
write.csv(schaefer200x7_networkpair_labels, "/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_netpair_labels_noduplicates.csv", row.names=F, quote=F)

num_networks = 17
schaefer200x17_networkpair <- read.csv(sprintf("%s/networkpair_subxnetpair_matrix_schaefer200x%s_orig.csv", measure_outputs_dir, num_networks))
schaefer200x17_networkpair_labels <- setdiff(names(schaefer200x17_networkpair), "subject")
write.csv(schaefer200x17_networkpair_labels, "/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_netpair_labels_noduplicates.csv", row.names=F, quote=F)


 

