library(dplyr)
library(gratia) 
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
library(NEST)

# a script to run NEST to see if sex differences are enriched for within-network connectivity of default mode network A and limbic temporal pole (from schaefer200x17)

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
metric = args[1] # i.e. within-network connectivity 
atlas = args[2]
network_of_interest = args[3] # i.e. DefaultA
dataset = "all_datasets"
print(paste("Running NEST for", dataset, metric))

################## 
# Set Directories 
################## 
outputs_root <- "/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets/"
NEST_outputs_dir <- paste0(outputs_root, "/NEST/", metric)
metric_outputs_dir <- paste0(outputs_root, "/", metric)

if (dir.exists(NEST_outputs_dir)) {
  print(paste(NEST_outputs_dir, "already exists"))
} else {
  dir.create(NEST_outputs_dir, recursive = TRUE)
  print(paste(NEST_outputs_dir, "created"))
}

################### 
# Define functions 
###################
# a wrapper function for running NEST  
NEST_wrapper <- function(metric, df, parcel.labels, network_of_interest) {
  # make sure demographics and data have same subject ordering
  if(!identical(as.character(demographics$sub), as.character(df$sub))) {
    stop("Error: The 'sub' columns in 'demographics' and 'df' do not match.")
  }
  df <- df %>% select(-sub)
  df <- df %>% select(-any_of("X")) # in case the dataframe has a column "X" (the covbat-ed datasets have this)

  # make matrix
  covs <- as.matrix(cbind( demographics$age, demographics$sex, as.numeric(demographics$meanFD_avgSes))) # matrix dimensions = N x 3
  colnames(covs) <- c("age", "sex", "meanFD_avgSes")
  
  # get column names
  colnames_df <- names(df)
  
  # final variables for NEST. See https://github.com/smweinst/NEST/blob/main/R/readme.md 
  X <- as.matrix(df)
  colnames(X) <- NULL
  dat <- covs

  net <- list(net_interest <- parcel.labels$network_of_interest) # list of regions that are in my network of interest
  cbind(names(df), unlist(net))
  
  print(paste("Running NEST for", metric))

  # one-sided: is delta adj rsq of sex more extreme inside vs. outside the network-of-interest?
  print("Running one-sided NEST")
  result_onesided <- NEST(
    statFun = "gam.deltaRsq",
    args = list(X = X, 
                dat = dat, 
                gam.full.formula = as.formula(X.v ~ s(age, k = 3, fx = TRUE) + sex + meanFD_avgSes), 
                gam.null.formula = as.formula(X.v ~ s(age, k = 3, fx = TRUE) + meanFD_avgSes), 
                lm.formula = as.formula(X.v ~ age + sex + meanFD_avgSes),  
                y.in.gam = "sex", 
                y.in.lm = "sex", 
                y.permute = "sex",  
                n.perm = 10000),
    net.maps = net,
    one.sided = TRUE,
    n.cores = 4, 
    seed = 123, 
    what.to.return = "everything"
  )
  saveRDS(result_onesided, paste0(NEST_outputs_dir, "/", "NEST_", network_of_interest, "_", metric, ".RData"))
  
}

################### 
# Load files  
###################
PNC_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "PNC"))
HCPD_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HCPD"))
NKI_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "NKI"))
HBN_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HBN"))

# load atlas labels for the relevant metric
if (metric %in% c("GBC", "BNC", "WNC")) {
  # load parcel labels
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/parcellations/%1$s_regionlist_final.csv", atlas))

} else if (metric == "edge") {
  # load edge labels 
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/edge/%1$s_edge.csv", atlas))
  parcel.labels <- parcel.labels[,2] 
  parcel.labels <- gsub("[0-9]+Networks", "Networks", parcel.labels)
  parcel.labels <- data.frame(parcel.labels)
  names(parcel.labels) <- "label"
  
} else { # network pairs
  # load network-pair labels
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/parcellations/%1$s_netpair_labels_noduplicates.csv", atlas))
  names(parcel.labels) <- "label"
}

# set network of interest
parcel.labels$network_of_interest <- ifelse(grepl(network_of_interest, parcel.labels$label), 1, 0)

# load demographics for all datasets and combine
PNC_demographics <- read.csv(PNC_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "PNC") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
HCPD_demographics <- read.csv(HCPD_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = paste0("HCPD_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
NKI_demographics <- read.csv(NKI_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "NKI") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
HBN_demographics <- read.csv(HBN_config_data$demographics, stringsAsFactors = TRUE) %>% rename(site = ses) %>% mutate(dataset_site = paste0("HBN_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
demographics <- rbind(PNC_demographics, HCPD_demographics, NKI_demographics, HBN_demographics)

# load functional connectivity data 
if(metric %in% c("GBC", "BNC", "WNC")) {
  filename <-  sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  data <- read.csv(sprintf("%1$s/%2$s_covbat.csv", metric_outputs_dir, filename))
} else if (metric == "edge") {
  data <- readRDS(sprintf("%s/subx%s_%s_covbat.RData", metric_outputs_dir, metric, atlas))
} else if (metric == "networkpair") {
  filename <-  sprintf("%s_%s_%s", metric, "subxnetpair_matrix", atlas) 
  data <- read.csv(sprintf("%1$s/%2$s_covbat.csv", metric_outputs_dir, filename))
}

data <- data %>% rename(sub = subject) 
names(data) <-  gsub("[A-Z0-9]+Networks", "Networks", names(data)) # remove leading numbers in colnames

################### 
# Run NEST
################### 
NEST_wrapper(metric, data, parcel.labels, network_of_interest)
