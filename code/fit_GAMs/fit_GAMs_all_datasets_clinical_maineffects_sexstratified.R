library(data.table)
library(dplyr)
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
source("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/fit_GAMs/gam_functions/GAM_functions_covariates.R")

# This script fits developmental regionwise GAMs on functional connectivity data using functions from GAM_functions_covariates.R 
# for the aggregated dataset that includes PNC, NKI, and HBN. 

# fyi: GAMs are fit for 1 specified metric  
# Specifically this script does the following:

# 1) run gam.fit.covariate to look at the main effects of clinical measures (interested in internalizing) on FC measures

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
metric = args[1]
atlas = args[2] # should be schaefer200x17
dataset = "all_datasets" # looking at the aggregated dataset

print(paste("Fitting sex_diff GAMs for", dataset, atlas))
print(paste("metric:", metric))

################## 
# Set directories  
################## 
outputs_root <- "/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets/"
GAM_outputs_dir <- paste0(outputs_root, "/GAM/", metric)
metric_outputs_dir <- paste0(outputs_root, "/", metric)

if (dir.exists(GAM_outputs_dir)) {
  print(paste(GAM_outputs_dir, "already exists"))
} else {
  dir.create(GAM_outputs_dir, recursive = TRUE)
  print(paste(GAM_outputs_dir, "created"))
}

PNC_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "PNC"))
HCPD_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HCPD"))
NKI_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "NKI"))
HBN_config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HBN"))

################### 
# Define functions 
###################
# load pfactor data and merge with demographics
load_pfactor_data <- function(file_path) {
  pfactor <- fread(file_path)
  pfactor <- pfactor %>% rename(sub = participant_id)
  pfactor$sub <- paste0("sub-", pfactor$sub)
  pfactor <- pfactor %>% select(sub, p_factor_mcelroy_harmonized_all_samples, internalizing_mcelroy_harmonized_all_samples, externalizing_mcelroy_harmonized_all_samples)
  pfactor <- pfactor[!duplicated(pfactor), ]
  pfactor <- pfactor[complete.cases(pfactor), ]
  return(pfactor)
}

# 1) run gam.fit.covariate to look at the main effect  
run_gam.fit.covariate <- function(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx, filename, atlas){
  sex.maineffects <- mclapply(parcel.labels, 
                              function(x){ as.data.frame(gam.fit.clinical.measures(gam.data = gam_df,  
                                                                          region = as.character(x), 
                                                                          smooth_var = smooth.var, 
                                                                          covariates.interest = covariates.interest, 
                                                                          covariates.noninterest = covariates.noninterest, 
                                                                          knots = k, 
                                                                          set_fx = set.fx)) %>% mutate(region = x)}, mc.cores = 4) 
  sex.maineffects <- do.call(rbind, sex.maineffects)
  write.csv(sex.maineffects, sprintf("%1$s/%2$s_GAM_%3$s_maineffects_%4$s.csv", GAM_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
}

 

################### 
# Load files  
###################
# load atlas labels for the relevant metric
if (metric %in% c("GBC", "BNC", "WNC")) {
  # load parcel labels
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/parcellations/%1$s_regionlist_final.csv", atlas))
  parcel.labels <- parcel.labels$label
  
} else if (metric == "edge") {
  # load edge labels 
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/edge/%1$s_edge.csv", atlas))
  parcel.labels <- parcel.labels[,2] 
  parcel.labels <- gsub("[0-9]+Networks", "Networks", parcel.labels)
  
} else { # network pairs
  # load network-pair labels
  parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/parcellations/%1$s_netpair_labels_noduplicates.csv", atlas))
  parcel.labels <- parcel.labels[,1]
}

# load demographics for all datasets and combine
PNC_demographics <- read.csv(PNC_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "PNC") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
NKI_demographics <- read.csv(NKI_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "NKI") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
HBN_demographics <- read.csv(HBN_config_data$demographics, stringsAsFactors = TRUE) %>% rename(site = ses) %>% mutate(dataset_site = paste0("HBN_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)

#load clinical measures
PNC_pfactor <- load_pfactor_data("/cbica/projects/network_replication/covariate_analyses/pfactor/input/PNC/sample_info/study-PNC_desc-participants.tsv")
NKI_pfactor <- load_pfactor_data("/cbica/projects/network_replication/covariate_analyses/pfactor/input/NKI/sample_info/study-NKI_desc-participants.tsv")
NKI_pfactor$sub <- gsub("sub-", "", NKI_pfactor$sub)
HBN_pfactor <- load_pfactor_data("/cbica/projects/network_replication/covariate_analyses/pfactor/input/HBN/sample_info/study-HBN_desc-participants.tsv")

PNC_demographics <- inner_join(PNC_pfactor, PNC_demographics,by = "sub")
NKI_demographics <- inner_join(NKI_pfactor, NKI_demographics,by = "sub")
HBN_demographics <- inner_join(HBN_pfactor, HBN_demographics,by = "sub")

demographics <- rbind(PNC_demographics, NKI_demographics, HBN_demographics)

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

# look at DMN only
#if(metric %in% c("GBC", "BNC", "WNC")) {
#  data <- data[ , grepl("Default", names(data)) | names(data) == "sub"]
#  parcel.labels <- parcel.labels[grepl("Default", parcel.labels)]
#} else if (metric == "edge") {
#  data <- data[, grepl("(Default.*Default)", names(data)) | names(data) == "sub"]
#  parcel.labels <- parcel.labels[grepl("(Default.*Default)", parcel.labels)]
  
#} else if (metric == "networkpair" ) {
#  data <- data[, grepl("(default.*default)", names(data)) | names(data) == "sub"]
#  parcel.labels <- parcel.labels[grepl("(default.*default)", parcel.labels)]
#}
  
  

#################################################### 
# Merge demographics and qc with FC data and sex_diff
#################################################### 
# merge
gam_df <- merge(data, demographics, by = "sub")
gam_df <- gam_df %>% drop_na(sex) 
write.table(gam_df$sub, paste0("/cbica/projects/network_replication/covariate_analyses/sex_diff/input/all_datasets//sample_info/final_sex_diff_clinical_sexstratified_sample.txt"), col.names=F, row.names=F, quote = TRUE)

smooth.var = "age"
k = 3
set.fx = T
covs = "meanFD_avgSes"
 
################### 
# Fit da GAMz 
###################
filenameF <- "clinical_sexstratifiedF"
filenameM <- "clinical_sexstratifiedM"

# 1) run gam.fit.covariate to look at the main effect of sex  
covariates.interest = c("p_factor_mcelroy_harmonized_all_samples", "internalizing_mcelroy_harmonized_all_samples","externalizing_mcelroy_harmonized_all_samples")
covariates.noninterest = "meanFD_avgSes"

gam_df_F <- gam_df %>% filter(sex == "Female")
gam_df_M <- gam_df %>% filter(sex == "Male")
run_gam.fit.covariate(gam_df_F, smooth.var, covariates.interest, covariates.noninterest, k, set.fx, filenameF, atlas)
run_gam.fit.covariate(gam_df_M, smooth.var, covariates.interest, covariates.noninterest, k, set.fx, filenameM, atlas)

print("Script finished!")
