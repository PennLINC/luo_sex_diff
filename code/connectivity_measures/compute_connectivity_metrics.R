library(cifti)
library(dplyr)
library(magrittr)
library(rjson)
library(stringr)

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
metric = args[2]
atlas = args[3] # i.e. "schaefer200x17"
num_networks = sub(".*x", "", atlas)

print(paste("Computing", metric, atlas))

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
outputs_root <- config_data$covariate_output_root
metric_output_dir <- paste0(outputs_root, "/", metric)

# outputs_root made in bash script
if (!dir.exists(metric_output_dir)) {
  # If directory doesn't exist, create it
  dir.create(metric_output_dir, recursive = TRUE)
  print(paste("Directory", metric_output_dir, "created."))
} else {
  print(paste("Directory", metric_output_dir, "already exists."))
}

################## 
# Define Functions
##################

# helper function for mapping parcels to network names  
map_network <- function(x, num_networks) {
  if (num_networks == 7) {
    dplyr::case_when(
      stringr::str_detect(x, "Vis")         ~ "Visual",
      stringr::str_detect(x, "SomMot")      ~ "SomMot",
      stringr::str_detect(x, "DorsAttn")    ~ "DorsAttn",
      stringr::str_detect(x, "SalVentAttn") ~ "SalVentAttn",
      stringr::str_detect(x, "Cont")        ~ "Frontoparietal",
      stringr::str_detect(x, "Limb")        ~ "Limbic",
      stringr::str_detect(x, "Default")     ~ "Default",
      TRUE                                  ~ x    # pass through original label if no match
    )
  } else if (num_networks == 17) {
    dplyr::case_when(
      stringr::str_detect(x, "VisCent")           ~ "visualCentral",
      stringr::str_detect(x, "VisPeri")           ~ "visualPeripheral",
      stringr::str_detect(x, "SomMotA")           ~ "somatomotorA",
      stringr::str_detect(x, "SomMotB")           ~ "somatomotorBAuditory",
      stringr::str_detect(x, "DorsAttnA")         ~ "dorsalAttentionA",
      stringr::str_detect(x, "DorsAttnB")         ~ "dorsalAttentionB",
      stringr::str_detect(x, "SalVentAttnA")      ~ "salienceVentralAttentionA",
      stringr::str_detect(x, "SalVentAttnB")      ~ "salienceVentralAttentionB",
      stringr::str_detect(x, "LimbicA_TempPole")  ~ "limbicTemporopolar",
      stringr::str_detect(x, "LimbicB_OFC")       ~ "limbicOrbitofrontal",
      stringr::str_detect(x, "ContA")             ~ "frontoparietalControlA",
      stringr::str_detect(x, "ContB")             ~ "frontoparietalControlB",
      stringr::str_detect(x, "ContC")             ~ "frontoparietalControlC",
      stringr::str_detect(x, "DefaultA")          ~ "defaultA",
      stringr::str_detect(x, "DefaultB")          ~ "defaultB",
      stringr::str_detect(x, "DefaultC")          ~ "defaultC",
      stringr::str_detect(x, "TempPar")           ~ "temporoparietal",
      TRUE                                  ~ x     # pass through original label if no match
    )
  }
} 


# Function for making empty dataframes for connectivity metric outputs
# @param subject_list list of participant ID's
# @param atlas A character string, name of atlas of interest (gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", or "schaefer400x17")
make_output_dfs <- function(subject_list, atlas){
  if(atlas == "glasser"){
    n <- 361
    regionheaders <- as.character(glasser.parcel.labels$label)
  } else if (atlas =="gordon"){
    n <- 334
    regionheaders <- as.character(gordon.parcel.labels$label)
  } else if (atlas =="schaefer200x7" | atlas=="schaefer200"){
    n <- 201
    regionheaders <- as.character(schaefer200x7.parcel.labels$label)
  } else if (atlas =="schaefer200x17") {
    n <- 201
    regionheaders <- as.character(schaefer200x17.parcel.labels$label)
  } else if (atlas =="schaefer400x7" | atlas=="schaefer400") {
    n <- 401
    regionheaders <- as.character(schaefer400x7.parcel.labels$label)
  } else if (atlas =="schaefer400x17") {
    n <- 401
    regionheaders <- as.character(schaefer400x17.parcel.labels$label)
  }
  subxparcel.matrix  <- matrix(data = NA, nrow = length(subject_list), ncol = n)
  demoheaders <- c("subject")
  colheaders <- as.matrix(c(demoheaders,regionheaders))
  colnames(subxparcel.matrix) <- colheaders
  return(subxparcel.matrix)
  
}

# functions for computing functional connectivity strength (FC_strength), 
# between-network coupling (BNC), within-network coupling (WNC), and edges

# Function for computing functional connectivity strength: average connectivity of a given region to every other region in the brain
# @param subject id of subject of interest
# @param atlas A character string, name of atlas of interest (e.g. glasser, gordon, schaefer200, or schaefer400)
# @param dataset A character string, name of dataset
computeFC_strength <- function(subject, atlas, dataset){
  
  #read in connectivity matrix
  if(dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") {
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData", dataset, subject))
    if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") {  
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]]
  } else {
    print("Provide valid dataset")
    
  }
   
  # compute average connectivity, excluding self-connectivity
  diag(connect.matrix) <- NA
  FC_strength <- as.array(rowMeans(connect.matrix, na.rm = TRUE))
  
  return(FC_strength)
}
 

# Function for Computing Between-Network or Within-Network Connectivity
# between-network connectivity: average connectivity of a given region to every other region in other networks
# within-network connectivity: average connectivity of a given region to every other region in its own network
# @param subject id of subject of interest
# @param atlas A character string, name of atlas of interest (gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", or "schaefer400x17")
# @param metric A character string, name of connectivity metric (either "BNC" or "WNC")
# @param dataset A character string, name of dataset
computeBNC_WNC <- function(subject, atlas, metric, dataset) {
  message(sprintf("Loading %s connectivity matrix...", atlas))
  
  # read in connectivity matrix
  #read in connectivity matrix 
  if (dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") { 
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData",dataset, subject)) 
    if (str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17") 
    } else { 
        atlas_name <- atlas 
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]] 
  } else if(dataset == "NKI") { 
      connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject)) 
      ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1] 
      if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) { 
        atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17") 
      } else { 
          atlas_name <- atlas 
      } 
      connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]] 
  } else { 
        print("Provide valid dataset") }
  
  # load parcel labels
  if (str_detect(atlas, "schaefer200x7")) {
    parcel.labels <- schaefer200x7.parcel.labels
    parcel.labels <- cbind(parcel.labels, map_network(parcel.labels$label, 7)) |>
      setNames(c("label", "network"))
  } else if (str_detect(atlas, "schaefer200x17")) {
    parcel.labels <- schaefer200x17.parcel.labels
    parcel.labels <- cbind(parcel.labels, map_network(parcel.labels$label, 17)) |>
      setNames(c("label", "network"))
  } else if (str_detect(atlas, "schaefer400x7")) {
    parcel.labels <- schaefer400x7.parcel.labels
    parcel.labels <- cbind(parcel.labels, map_network(parcel.labels$label, 7)) |>
      setNames(c("label", "network"))
  } else if (str_detect(atlas, "schaefer400x17")) {
    parcel.labels <- schaefer400x17.parcel.labels
    parcel.labels <- cbind(parcel.labels, map_network(parcel.labels$label, 17)) |>
      setNames(c("label", "network"))
  } else {
    stop("Please provide a valid atlas: schaefer200x7, schaefer200x17, schaefer400x7, or schaefer400x17.")
  }
  
  # make sure matrix is right size
  stopifnot(nrow(connect.matrix) == ncol(connect.matrix),
            nrow(connect.matrix) == nrow(parcel.labels))
  
  # ensure symmetry and remove diagonals  
  connect.sym <- (connect.matrix + t(connect.matrix)) / 2
  diag(connect.sym) <- NA
  
  networks <- parcel.labels$network
  
  # compute metric
  if (metric == "BNC") {
    mask <- outer(networks, networks, FUN = "!=")
  } else if (metric == "WNC") {
    mask <- outer(networks, networks, FUN = "==")
    diag(mask) <- FALSE
  } else {
    stop("metric must be 'BNC' or 'WNC'.")
  }
  
  mat.masked <- connect.sym
  mat.masked[!mask] <- NA
  FC_metric <- rowMeans(mat.masked, na.rm = TRUE)
  
  return(FC_metric)
}

  
# Function for Extracting Parcel-Parcel Connectivity (edges)
# @param subject rbcid of subject of interest
# @param atlas A character string, name of atlas of interest (e.g. glasser, gordon, schaefer200x7, schaefer200x17, schaefer400x7, or schaefer400x17)
# @param dataset A character string, name of dataset
extractParcel2ParcelConn <- function(subject, atlas, dataset){
  
  #read in connectivity matrix  
  if(dataset == "PNC" | dataset=="HCPD" | dataset=="HBN") {
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData", dataset, subject))
    if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") { #check if BAS1 exists for subject
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(str_detect(atlas, "schaefer")) {
      atlas_name <- str_extract(atlas, "schaefer[0-9]")
      atlas_name <- paste0(atlas_name, "17")
    } else {
      atlas_name <- atlas
    }
    matrixName <- paste0(ses_name, "_", atlas_name, "_conn")
    connect.matrix <- connect.matrix[[matrixName]]
  }  
    
  if (atlas == "glasser"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/glasser_edge.csv")
  } else if (atlas == "gordon"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/gordon_edge.csv")
  }  else if(atlas == "schaefer200x7"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/schaefer200x7_edge.csv")
  } else if(atlas == "schaefer400x7"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/schaefer400x7_edge.csv")
  } else if(atlas == "schaefer200x17"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/schaefer200x17_edge.csv")
  } else if(atlas == "schaefer400x17"){
    edge <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/schaefer400x17_edge.csv")
  } else {
    print("Please provide valid atlas (glasser, gordon, schaefer200, or schaefer400")
  } 
   
  # set upper triangle of matrix to 9999 and vectorize the matrix
  mat <- connect.matrix
  mat[upper.tri(mat)] <- 9999
  connect.vector <- c(mat)
  connect.vector <- connect.vector[-c(which(connect.vector==9999))] # this doesn't remove the diagonal
  connect.vector <- connect.vector[-c(which(connect.vector==1))] # remove edges that connect the same parcel to itself 
  print(paste(which(participants %in% subject), "/", length(participants), "-", atlas))
  return(connect.vector)
  
}


# Helper function for computing average connectivity between each pair of networks and within networks
compute_mean_networkpair <- function(
    mat,
    map_network,
    net_order    = NULL,
    exclude_diag = FALSE,   # NOTE: default changed to match collapse_to_network()
    na_rm        = TRUE
) {
  stopifnot(is.matrix(mat),
            !is.null(rownames(mat)), !is.null(colnames(mat)))
  
  # map parcels -> networks using name patterns
  r_comm <- map_network(rownames(mat), num_networks)
  c_comm <- map_network(colnames(mat), num_networks)
  
  # safety: if mapper ever returned NA, fall back to the original label (like TRUE ~ x does)
  if (anyNA(r_comm)) r_comm[is.na(r_comm)] <- rownames(mat)[is.na(r_comm)]
  if (anyNA(c_comm)) c_comm[is.na(c_comm)] <- colnames(mat)[is.na(c_comm)]
  
  # indices per network
  row_idx <- split(seq_len(nrow(mat)), r_comm)
  col_idx <- split(seq_len(ncol(mat)), c_comm)
  
  # choose order (union of present row/col network names)
  if (is.null(net_order)) {
    net_order <- union(names(row_idx), names(col_idx))
    net_order <- unique(net_order)
  }
  row_idx <- row_idx[net_order]
  col_idx <- col_idx[net_order]
  
  # compute block means
  net_mat <- matrix(NA_real_,
                    nrow = length(row_idx), ncol = length(col_idx),
                    dimnames = list(names(row_idx), names(col_idx)))
  
  for (i in seq_along(row_idx)) {
    for (j in seq_along(col_idx)) {
      block <- mat[row_idx[[i]], col_idx[[j]], drop = FALSE]
      # optional: drop within-network parcel self-terms if requested
      if (exclude_diag && names(row_idx)[i] == names(col_idx)[j]) {
        diag(block) <- NA_real_
      }
      net_mat[i, j] <- mean(block, na.rm = na_rm)
    }
  }
  net_mat
}

# Function for formatting output for network pairs
# One row per subject; columns for every unique (i <= j) network pair
net_pairs_wide_all <- function(mat, subject, sep = "_", net_order = NULL) {
  stopifnot(is.matrix(mat),
            !is.null(rownames(mat)), !is.null(colnames(mat)))
  
  # Optionally reorder to a desired network order for stable column layout
  if (!is.null(net_order)) {
    miss_r <- setdiff(net_order, rownames(mat))
    miss_c <- setdiff(net_order, colnames(mat))
    if (length(miss_r) || length(miss_c)) {
      stop("net_order labels not found in matrix dimnames: ",
           paste(unique(c(miss_r, miss_c)), collapse = ", "))
    }
    mat <- mat[net_order, net_order, drop = FALSE]
  }
  
  # Take upper triangle (incl. diagonal) to keep each unordered pair once
  inds <- upper.tri(mat, diag = TRUE)
  ij   <- which(inds, arr.ind = TRUE)
  rn   <- rownames(mat)[ij[, 1]]
  cn   <- colnames(mat)[ij[, 2]]
  keys <- paste(rn, cn, sep = sep)
  vals <- mat[inds]
  
  # Build a one-row data frame: subject + pair columns
  out <- as.data.frame(t(vals), stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- keys
  tibble(subject = subject) |> bind_cols(out)
}


# Function for computing average connectivity between and within networks for each network pair
# @param subject String for id of subject of interest
# @param atlas String for name of atlas 
wrapper_compute_mean_networkpair <- function(subject, atlas) {
  print(subject)
  parcel.labels <- get(paste0(atlas, ".parcel.labels"))
  
  # read in connectivity matrix (unchanged from your code)
  if (dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") {
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData", dataset, subject))
    
    if (str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(stringr::str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    }
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
    
  } else if (dataset == "NKI") {
    connect.matrix <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/manuscript/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- stringr::str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if (str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(stringr::str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    }
    connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]]
    
  } else {
    stop("Provide valid dataset")
  }
  
  # ensure parcel names are present
  colnames(connect.matrix) <- parcel.labels$label
  rownames(connect.matrix) <- parcel.labels$label
  
  # >>> key change: use the mapper-based aggregator <<<
  mean_network_conn <- compute_mean_networkpair(
    connect.matrix,
    map_network   = map_network,
    exclude_diag  = FALSE,  # match collapse_to_network default
    na_rm         = TRUE
  )
  
  mean_network_conn_wide <- net_pairs_wide_all(mean_network_conn, subject)
  print(paste(which(participants %in% subject), "/", length(participants), "-", subject))
  return(mean_network_conn_wide)
}

################## 
# Read files 
################## 
# load atlases (for extractParcel2ParcelConn)
glasser.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")
gordon.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/gordon_regionlist_final.csv")
schaefer200x7.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
schaefer200x17.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv") # the atlas used in this project
schaefer400x7.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")
schaefer400x17.parcel.labels <- read.csv("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")

# load participant list
participants <- read.table(config_data$subject_list)[[1]]
if (dataset=="NKI") {
  participants <- gsub("sub-", "", participants)
} else {
  participants <- participants
}

################## 
# Compute metric 
################## 

if (metric == "FC_strength"){
  # make empty dataframes for connectivity metric outputs
  subxparcel.matrix <- make_output_dfs(participants, atlas)
  
  # compute FC_strength
  for(sub in c(1:length(participants))){
    subjectID=as.character(participants[sub])
    print(paste(subjectID, atlas, metric, dataset))
    id.data <- computeFC_strength(subjectID, atlas, dataset)
    df_toUpdate <- subxparcel.matrix
    df_toUpdate[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to subxparcel.matrix
    assign("subxparcel.matrix", df_toUpdate)  
    print(paste(sub, "/", length(participants), "-", subjectID, atlas))
  }
  
  # save out
  filename <- sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  write.csv(subxparcel.matrix, sprintf("%1$s/%2$s.csv", metric_output_dir, filename), row.names=F, quote=F)
  
} else if(metric == "BNC" | metric == "WNC") {
  # make empty dataframes for connectivity metric outputs
  subxparcel.matrix <- make_output_dfs(participants, atlas)
  
  # compute BNC or WNC
  for(sub in c(1:length(participants))){
    subjectID=as.character(participants[sub])
    print(paste(subjectID, atlas, metric, dataset))
    id.data <- computeBNC_WNC(subjectID, atlas, metric, dataset)
    
    df_toUpdate <- subxparcel.matrix
    df_toUpdate[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to subxparcel.matrix
    assign("subxparcel.matrix", df_toUpdate)  
    print(paste(sub, "/", length(participants), "-", subjectID, atlas))
  }
  
  # save out
  filename <- sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  write.csv(subxparcel.matrix, sprintf("%1$s/%2$s.csv", metric_output_dir, filename), row.names=F, quote=F)
  
} else if(metric == "edge") {
  # extract parcel-parcel connectivity
  extractEdge <- lapply(participants, extractParcel2ParcelConn, atlas, dataset)
  names(extractEdge) <- participants 
  columns_extractEdge <- bind_rows(extractEdge) # turn list into a dataframe
  subxedge <- as.data.frame(t(columns_extractEdge)) # transpose dataframe
  
  # load edge labels
  edge <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/%1$s_edge.csv", atlas))
  edge <- edge[,2]
  names(subxedge) <- edge # columns = names of edges
  
  # format subxedge df
  subxedge <- subxedge %>% mutate(subject = rownames(subxedge))
  subxedge <- subxedge %>% relocate(subject)
  rownames(subxedge) <- NULL
  
  # save out
  saveRDS(subxedge, sprintf("%1$s/subxedge_%2$s.RData", metric_output_dir, atlas)) 
  print(paste("Edges for", atlas, "done"))
  
} else if(metric == "networkpair") {
  subjects_mean_netconn <- lapply(participants, wrapper_compute_mean_networkpair,  atlas=atlas)
  names(subjects_mean_netconn) <- participants
  
  saveRDS(subjects_mean_netconn, sprintf("%s/networkpair_subxnetpair_matrix_%s_orig_notbound.RData", metric_output_dir, atlas))
  subjects_mean_netconn <- bind_rows(subjects_mean_netconn)
  
  # creates a subject x network_pair csv (nrow=length(subjects), ncol=length(network_pair))
  write.csv(subjects_mean_netconn, sprintf("%s/networkpair_subxnetpair_matrix_%s_orig.csv", metric_output_dir, atlas), row.names=F, quote=F)
  
}
 