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
dataset = args[1]
metric = args[2]
atlas = args[3] # i.e. schaefer200x7
print(paste("Processing", dataset, metric))

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
outputs_root <- config_data$covariate_output_root
metric_outputs_dir <- paste0(outputs_root, "/", metric)

# load demographics
demographics <- read.csv(config_data$demographics, stringsAsFactors = TRUE)
if (dataset == "HBN") {
  demographics <- demographics %>% rename(site = ses)
}

################## 
# Load files 
################## 
if(metric %in% c("FC_strength", "BNC", "WNC")) {
  filename <-  sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  subxparcel <- read.csv(sprintf("%1$s/%2$s.csv", metric_outputs_dir, filename))
} else if (metric == "edge") {
  subxparcel <- readRDS(sprintf("%s/subx%s_%s.RData", metric_outputs_dir, metric, atlas))
} else if (metric == "networkpair") {
  subxparcel <- read.csv(sprintf("%s/networkpair_subxnetpair_matrix_%s_orig.csv", metric_outputs_dir, atlas))
}

print(paste(metric, "loaded"))

subxparcel <- subxparcel[complete.cases(subxparcel),]
names(subxparcel) <- gsub("\\.", "-", gsub("^X", "", names(subxparcel)))

if (identical(as.factor(demographics$sub), as.factor(subxparcel$subject))) {
  subxparcel <- data.frame(subxparcel)[,-1]
} else {
  stop("ERROR: Subject IDs in demographics and subxparcel do not match. Exiting.")
}

row.names(subxparcel) <- demographics$sub

# set covariate vectors
age_vec <- demographics$age 
sex_vec <- as.factor(demographics$sex) 
meanFD_avgSes_vec <- demographics$meanFD_avgSes 
covar_df <- bind_cols(demographics$sub, as.numeric(age_vec), as.factor(sex_vec), as.numeric(meanFD_avgSes_vec))
covar_df <- dplyr::rename(covar_df, sub=...1,
                          age = ...2,
                          sex = ...3,
                          meanFD_avgSes = ...4)

################## 
# Harmonize data 
################## 
data.harmonized <- covfam(subxparcel, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(meanFD_avgSes))
print(paste0(metric, " harmonized"))

# clean up covbat output for saving to RData
subxparcel_covbat <- data.frame(data.harmonized$dat.covbat)
subxparcel_covbat$subject <- rownames(subxparcel)
subxparcel_covbat <- subxparcel_covbat %>% relocate(subject)
rownames(subxparcel_covbat) <- NULL

# save out!
if(metric != "edge") {
  filename <-  sprintf("%s_%s_%s", metric, "subxparcel_matrix", atlas) 
  write.csv(subxparcel_covbat, sprintf("%1$s/%2$s_covbat.csv", metric_outputs_dir, filename))
} else if (metric == "edge") {
  saveRDS(subxparcel_covbat, sprintf("%s/subx%s_%s_covbat.RData", metric_outputs_dir, metric, atlas))
}

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

print("File saved")


 