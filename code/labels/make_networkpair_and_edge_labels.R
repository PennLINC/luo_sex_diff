library(rjson)


# Set project root depending on whether code is running on PARCC or locally
if (dir.exists("/ceph/projects/sattertt/pennlinc-parcc/network_replication")) {
  root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication"
} else if (dir.exists("~/parcc/network_replication")) {
  root <- path.expand("~/parcc/network_replication")
} else {
  stop("Could not find network_replication project.")
}


################################################## 
# ------------ network pair labels --------------#
################################################## 

# Set Variables 
dataset = "PNC" # any dataset works
measure = "networkpair"
 
# Set Directories 
config_data <- fromJSON(file=sprintf("%1$s/covariate_analyses/sex_diff/code/config/config_%2$s.json", root, dataset))
outputs_root <- file.path(sprintf("%1$s/covariate_analyses/sex_diff/output/%2$s", root, dataset))
measure_outputs_dir <- paste0(outputs_root, "/", measure)
 
# Write out labels for each pair of networks
num_networks = 17
schaefer200x17_networkpair <- read.csv(sprintf("%s/networkpair_subxnetpair_matrix_schaefer200x%s_orig.csv", measure_outputs_dir, num_networks))
schaefer200x17_networkpair_labels <- setdiff(names(schaefer200x17_networkpair), "subject")
write.csv(schaefer200x17_networkpair_labels, sprintf("%1$s/atlases/parcellations/schaefer200x17_netpair_labels_noduplicates.csv", root), row.names=F, quote=F)



########################################## 
# ------------ edge labels --------------#
########################################## 

# function for creating vector of edge names, excluding self-connecting edges
# @param nrows_var Number of rows in parcel_labels
# @param atlas_edges An empty vector to be populated by edge names
# @param parcel_labels A df of parcel labels
make_edges <- function(nrows_var, atlas_edges, parcel_labels) {
  for (i in seq_len(nrows_var - 1)) {
    for (j in (i + 1):nrows_var) {
      x <- paste0(parcel_labels$orig_parcelname[i], "_to_", parcel_labels$orig_parcelname[j])
      atlas_edges <- append(atlas_edges, x)
      print(paste(i, j))
    }
  }
  return(atlas_edges)
}

schaefer200x17.parcel.labels <- read.csv(sprintf("%1$s/atlases/parcellations/schaefer200x17_regionlist_final.csv", root))
names(schaefer200x17.parcel.labels)[1] <- "orig_parcelname"

rows_schaefer200x17 <- nrow(schaefer200x17.parcel.labels)
schaefer200x17_edges <- c()
schaefer200x17_edges <- make_edges(rows_schaefer200x17, schaefer200x17_edges, schaefer200x17.parcel.labels)
schaefer200x17_edges <- gsub("Network", "17Network", schaefer200x17_edges)
write.csv(schaefer200x17_edges, sprintf("%1$s/atlases/edge/schaefer200x17_edge.csv", root), row.names=F, quote=F)