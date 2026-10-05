library(ComBatFamily)
library(data.table)
library(dplyr)
library(mgcv)
library(rjson)
library(stringr)
library(tidyr)

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
metric = args[1]
atlas = args[2] # i.e. schaefer200x17
print(paste("Processing", metric, atlas))
print(paste("Harmonizing across datasets for", metric, atlas))

 
################## 
# Set Directories 
################## 
PNC_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "PNC"))
HCPD_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HCPD"))
NKI_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "NKI"))
HBN_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HBN"))


outputs_root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/output/all_datasets"
metric_outputs_dir <- paste0(outputs_root, "/", metric)
PNC_outputs_dir <- paste0(PNC_config_data$covariate_output_root, "/", metric)
HCPD_outputs_dir <- paste0(HCPD_config_data$covariate_output_root, "/", metric)
NKI_outputs_dir <- paste0(NKI_config_data$covariate_output_root, "/", metric)
HBN_outputs_dir <- paste0(HBN_config_data$covariate_output_root, "/", metric)


if (!dir.exists(metric_outputs_dir)) {
  # If directory doesn't exist, create it
  dir.create(metric_outputs_dir, recursive = TRUE)
  print(paste("Directory", metric_outputs_dir, "created."))
} else {
  print(paste("Directory",metric_outputs_dir, "already exists."))
}


################## 
# Load files 
################## 
# load demographics for all datasets
PNC_demographics <- read.csv(PNC_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "PNC") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
HCPD_demographics <- read.csv(HCPD_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = paste0("HCPD_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
NKI_demographics <- read.csv(NKI_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "NKI") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
HBN_demographics <- read.csv(HBN_config_data$demographics, stringsAsFactors = TRUE) %>% rename(site = ses) %>% mutate(dataset_site = paste0("HBN_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
 
# combine demographics
all_datasets_demographics <- rbind(PNC_demographics, HCPD_demographics, NKI_demographics, HBN_demographics)

# load subxparcel files for all datasets  
if(metric %in% c("FC_strength", "BNC", "WNC")) {
  filename <-  sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  PNC_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", PNC_outputs_dir, filename))
  HCPD_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", HCPD_outputs_dir, filename))
  NKI_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", NKI_outputs_dir, filename))
  HBN_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", HBN_outputs_dir, filename))
} else if (metric == "edge") {
  PNC_subxparcel <- readRDS(sprintf("%s/subx%s_%s.RData", PNC_outputs_dir, metric, atlas))
  HCPD_subxparcel <- readRDS(sprintf("%s/subx%s_%s.RData", HCPD_outputs_dir, metric, atlas))
  NKI_subxparcel <- readRDS(sprintf("%s/subx%s_%s.RData", NKI_outputs_dir, metric, atlas))
  HBN_subxparcel <- readRDS(sprintf("%s/subx%s_%s.RData", HBN_outputs_dir, metric, atlas))
} else if (metric == "networkpair") {
  filename <-  sprintf("networkpair_subxnetpair_matrix_%s_orig", atlas) 
  PNC_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", PNC_outputs_dir, filename))
  HCPD_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", HCPD_outputs_dir, filename))
  NKI_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", NKI_outputs_dir, filename))
  HBN_subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", HBN_outputs_dir, filename))
}

print(paste(metric, "loaded"))

all_datasets_subxparcel <- rbind(PNC_subxparcel, HCPD_subxparcel, NKI_subxparcel, HBN_subxparcel)
all_datasets_subxparcel <- all_datasets_subxparcel[complete.cases(all_datasets_subxparcel),]
names(all_datasets_subxparcel) <- gsub("\\.", "-", gsub("^X", "", names(all_datasets_subxparcel)))

if (identical(as.character(all_datasets_demographics$sub), as.character(all_datasets_subxparcel$subject))) {
  all_datasets_subxparcel <- data.frame(all_datasets_subxparcel)[,-1]
} else {
  stop("ERROR: Subject IDs in all_datasets_demographics and all_datasets_subxparcel do not match. Exiting.")
}

row.names(all_datasets_subxparcel) <- all_datasets_demographics$sub

# set covariate vectors
age_vec <- all_datasets_demographics$age 
sex_vec <- as.factor(all_datasets_demographics$sex) 
meanFD_avgSes_vec <- all_datasets_demographics$meanFD_avgSes 
covar_df <- bind_cols(all_datasets_demographics$sub, as.numeric(age_vec), as.factor(sex_vec), as.numeric(meanFD_avgSes_vec))
covar_df <- dplyr::rename(covar_df, sub=...1,
                          age = ...2,
                          sex = ...3,
                          meanFD_avgSes = ...4)

################## 
# Harmonize data 
################## 
data.harmonized <- covfam(all_datasets_subxparcel, bat = as.factor(all_datasets_demographics$dataset_site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(meanFD_avgSes))
print(paste0(metric, "harmonized"))

# clean up covbat output for saving to RData
subxparcel_covbat <- data.frame(data.harmonized$dat.covbat)
subxparcel_covbat$subject <- rownames(all_datasets_subxparcel)
subxparcel_covbat <- subxparcel_covbat %>% relocate(subject)
rownames(subxparcel_covbat) <- NULL

# save out!
if(metric %in% c("FC_strength", "BNC", "WNC")) {
  filename <-  sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  write.csv(subxparcel_covbat, sprintf("%1$s/%2$s_covbat.csv", metric_outputs_dir, filename))
} else if (metric == "edge") {
  saveRDS(subxparcel_covbat, sprintf("%s/subx%s_%s_covbat.RData", metric_outputs_dir, metric, atlas))
} else if (metric == "networkpair") {
  filename <-  sprintf("%s_%s_%s", metric, "subxnetpair_matrix", atlas) 
  write.csv(subxparcel_covbat, sprintf("%1$s/%2$s_covbat.csv", metric_outputs_dir, filename))
}
 