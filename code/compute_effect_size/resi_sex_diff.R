library(data.table)
library(dplyr)
library(mgcv)
library(parallel)
library(rjson)
library(RESI)
library(stringr)
library(tidyr)
source("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/compute_effect_size/resiPEse_Th1_functions.R")


# This script fits developmental edge-, region-, and network-level GAMs on functional connectivity data 
# fyi: GAMs are fit for 1 specified metric (i.e. FC_strength)
# then it computes an effect size using RESI (with the resiPEse_Th1 and t2S methods)!
 
################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
metric = args[2]
atlas = args[3]

print(paste("Fitting sex_diff GAMs and computing RESI for", dataset, atlas))
print(paste("metric:", metric))

################## 
# Set directories  
################## 
if (dataset != "all_datasets") {
  config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset)) 
  outputs_root <- config_data$covariate_output_root
} else if (dataset == "all_datasets") {
  PNC_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "PNC"))
  HCPD_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HCPD"))
  NKI_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "NKI"))
  HBN_config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", "HBN"))
  outputs_root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/output/all_datasets"
  
}

resi_outputs_dir <- paste0(outputs_root, "/resi/", metric)
metric_outputs_dir <- paste0(outputs_root, "/", metric)

if (dir.exists(resi_outputs_dir)) {
  print(paste(resi_outputs_dir, "already exists"))
} else {
  dir.create(resi_outputs_dir, recursive = TRUE)
  print(paste(resi_outputs_dir, "created"))
}

################### 
# Define functions
###################
# Use the theory-based resiPEse_Th1 method (fast CI, no bootstrap)
# Fit GAM (1 GAM per region or network-pair) and compute the sex effect size with resiPEse_Th1
fit_and_resiPEse_Th1 <- function(region) {
  # Extract the connectivity value (response variable) for this region (or network-pair)
  # Each column in gam_df corresponds to one region’s connectivity measure
  y <- gam_df[[region]]
  # Quality control: skip parcels with no usable data
  # - all NA values
  # - near-zero variance (model would be unstable / meaningless)
  if (all(is.na(y)) || sd(y, na.rm=TRUE) < 1e-8) {
    return(data.frame(region=region, RESI_Xinyu=NA, LCI=NA, UCI=NA))
  }
  # GAM formula:
  formula <- as.formula(paste0("y ~ s(age, k = ", k, ", fx = ", set.fx, ") + ", covariates.noninterest, " + ",covariate.interest))
  
  # fit GAM
  mod <- mgcv::gam(formula, data = gam_df, method = "REML")
  
  
  # Compute standardized effect size for the sex coefficient using RESI:
  # - variable="sexFemale" gives the signed female-vs-male contrast
  # - unsigned=FALSE preserves direction of the effect
  # - (positive RESI = higher values in females relative to the male reference level)
  # - type="HC0" uses a basic robust covariance estimator
  # - returns RESI estimate + theory-based SE and CI (no bootstrap)
  out <- resiPEse_Th1(mod, variable="sexFemale", unsigned=FALSE, torz="t", type="HC0")
  
  # return a one-row summary for this region (or network-pair):
  # region ID + effect size estimate + CI
  cbind(region=region, out)
}


# t2S method (point-estimate only, very fast, t2S() uses the mgcv t-statistic directly):
# Fit GAM (1 GAM per region or network-pair) and compute the sex effect size with RESI (t2S)
fit_and_t2S <- function(region) {
  y <- gam_df[[region]]
  
  # Quality control: skip features with no usable data
  # - all NA values
  # - near-zero variance (model would be unstable / meaningless)
  if (all(is.na(y)) || sd(y, na.rm = TRUE) < 1e-8) {
    return(data.frame(region = region, RESI = NA_real_, t = NA_real_, p = NA_real_))
  }
  
  # GAM formula
  formula <- as.formula(paste0("y ~ s(age, k = ", k, ", fx = ", set.fx, ") + ", covariates.noninterest, " + ",covariate.interest))
  
  # Fit GAM
  mod <- mgcv::gam(formula, data = gam_df, method = "REML")
  
  # Compute RESI point estimate for the sex effect
  # - uses the t-statistic from the GAM under the hood
  # - fast and stable for large-scale (e.g., edge-level) analyses
  resi_sex <- RESI::t2S(
    summary(mod)$p.table["sexFemale", "t value"],
    rdf = mod$df.residual,
    n   = nrow(gam_df),
    unbiased = TRUE
  )
  
  # Return a one-row summary for this region (or edge/network-pair):
  # region ID + RESI point estimate
  data.frame(region = region, RESI = as.numeric(resi_sex))
}


################### 
# Load files  
###################
# load atlas labels for the relevant metric
if (metric %in% c("FC_strength", "BNC", "WNC")) {
  # load parcel labels
  parcel.labels <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/%1$s_regionlist_final.csv", atlas))
  parcel.labels <- parcel.labels$label
  
} else if (metric == "edge") {
  # load edge labels 
  parcel.labels <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/%1$s_edge.csv", atlas))
  parcel.labels <- parcel.labels[,1] 
  parcel.labels <- gsub("[0-9]+Networks", "Networks", parcel.labels)
  
} else { # network pairs
  # load network-pair labels
  parcel.labels <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/%1$s_netpair_labels_noduplicates.csv", atlas))
  parcel.labels <- parcel.labels[,1]
}

# load demographics
if (dataset != "all_datasets") {
  demographics <- read.csv(config_data$demographics, stringsAsFactors = TRUE)
  demographics <- demographics %>% select(sub, sex, age, meanFD_avgSes)
  demographics$sex <- factor(demographics$sex, levels = c("Male", "Female"))
} else if (dataset == "all_datasets") {
  # load demographics for all datasets and combine
  PNC_demographics <- read.csv(PNC_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "PNC") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
  HCPD_demographics <- read.csv(HCPD_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = paste0("HCPD_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
  NKI_demographics <- read.csv(NKI_config_data$demographics, stringsAsFactors = TRUE) %>% mutate(dataset_site = "NKI") %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
  HBN_demographics <- read.csv(HBN_config_data$demographics, stringsAsFactors = TRUE) %>% rename(site = ses) %>% mutate(dataset_site = paste0("HBN_", site)) %>% select(sub, age, sex, meanFD_avgSes, dataset_site)
  demographics <- rbind(PNC_demographics, HCPD_demographics, NKI_demographics, HBN_demographics)
  demographics$sex <- factor(demographics$sex, levels = c("Male", "Female"))
}

 
# load functional connectivity data (accommodating for different file names for datasets that got harmonized)
covbat_datasets <- c("HCPD", "HBN", "all_datasets")
csv_datasets <- c("PNC", "NKI")

if (metric %in% c("FC_strength", "BNC", "WNC")) {
  filename <- sprintf("%s_subxparcel_matrix_%s", metric, atlas)
  filepath <- sprintf("%s/%s%s.csv", metric_outputs_dir, filename, if (dataset %in% covbat_datasets) "_covbat" else "")
  data <- read.csv(filepath)
  
} else if (metric == "networkpair") {
  filename <- sprintf("%s_subxnetpair_matrix_%s", metric, atlas)
  filepath <- sprintf("%s/%s%s.csv", metric_outputs_dir, filename, if (dataset %in% covbat_datasets) "_covbat" else "_orig")
  data <- read.csv(filepath)
  
} else if (metric == "edge") {
  suffix <- if (dataset %in% covbat_datasets) "_covbat" else ""
  filepath <- sprintf("%s/subx%s_%s%s.RData", metric_outputs_dir, metric, atlas, suffix)
  data <- readRDS(filepath)
}

data <- data %>% rename(sub = subject) 
names(data) <-  gsub("[A-Z0-9]+Networks", "Networks", names(data)) # remove leading numbers in colnames


#################################################### 
# Merge demographics and qc with FC data and sex_diff
#################################################### 
# merge
gam_df <- merge(data, demographics, by = "sub")
gam_df <- gam_df %>% drop_na(sex) 

# check for sex diff in meanFD_avgSes (average average mean_fd across concatenated sessions)
print(t.test(gam_df$meanFD_avgSes[gam_df$sex=="Male"], gam_df$meanFD_avgSes[gam_df$sex=="Female"]))

 
################### 
# Run RESI
###################
# set filename variable for saving out GAM results
filename = "sex"
k = 3
set.fx = T
covariate.interest = "sex"
covariates.noninterest = "meanFD_avgSes"
 

if (metric %in% c("FC_strength", "BNC", "WNC", "networkpair")) { # run both resiPEse_Th1 and t2S for comparison
  # run resiPEse_Th1
  print("running resiPEse_Th1")
  resiPEse_Th1_output <- parallel::mclapply(parcel.labels, function(r) {
    tryCatch(fit_and_resiPEse_Th1(r),
             error=function(e) data.frame(region=r, error=conditionMessage(e)))
  }, mc.cores=2)
  resiPEse_Th1_output <- do.call(rbind, resiPEse_Th1_output)
  write.csv(resiPEse_Th1_output, sprintf("%1$s/%2$s_resiPEse_Th1_%3$s_maineffects_%4$s.csv", resi_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
  
  
  # run t2S
  print("running t2S")
  t2S_output <- parallel::mclapply(parcel.labels, function(r) {
    tryCatch(fit_and_t2S(r),
             error=function(e) data.frame(region=r, error=conditionMessage(e)))
  }, mc.cores=2)
  t2S_output <- do.call(rbind, t2S_output)
  write.csv(t2S_output, sprintf("%1$s/%2$s_t2S_%3$s_maineffects_%4$s.csv", resi_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
  
} else if (metric == "edge") { # run only t2S (less computationally heavy for so many edges!)
  # run t2S
  print("running t2S")
  t2S_output <- parallel::mclapply(parcel.labels, function(r) {
    tryCatch(fit_and_t2S(r),
             error=function(e) data.frame(region=r, error=conditionMessage(e)))
  }, mc.cores=2)
  t2S_output <- do.call(rbind, t2S_output) # only outputs the RESI effect size
  write.csv(t2S_output, sprintf("%1$s/%2$s_t2S_%3$s_maineffects_%4$s.csv", resi_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
  
}
 
print("Script finished!")
 