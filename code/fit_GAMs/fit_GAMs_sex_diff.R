library(data.table)
library(dplyr)
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
source("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/fit_GAMs/gam_functions/GAM_functions_covariates.R")

# This script fits developmental regionwise GAMs on functional connectivity data using functions from GAM_functions_covariates.R for each dataset.
# fyi: GAMs are fit for 1 specified metric (i.e. FC_strength)
# Specifically this script does the following:

# 1) run gam.fit.covariate to look at the main effect of sex  
# 2) run gam.factorsmooth.interaction to look at the sex-by-age interaction 
# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectory for female and males

# adapted from https://github.com/PennLINC/thalamocortical_development/blob/main/gam_functions/fit_envGams.R
# and https://github.com/PennLINC/luo_wm_dev/blob/main/code/fit_GAMs/fit_GAMs_development.R 

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
metric = args[2]
atlas = args[3]

print(paste("Fitting sex_diff GAMs for", dataset, atlas))
print(paste("metric:", metric))

################## 
# Set directories  
################## 
config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
outputs_root <- config_data$covariate_output_root
GAM_outputs_dir <- paste0(outputs_root, "/GAM/", metric)
metric_outputs_dir <- paste0(outputs_root, "/", metric)

if (dir.exists(GAM_outputs_dir)) {
  print(paste(GAM_outputs_dir, "already exists"))
} else {
  dir.create(GAM_outputs_dir, recursive = TRUE)
  print(paste(GAM_outputs_dir, "created"))
}

################### 
# Define functions 
###################
# 1) run gam.fit.covariate to look at the main effect of sex  
run_gam.fit.covariate <- function(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx, filename, atlas){
  sex.maineffects <- mclapply(parcel.labels, 
                              function(x){as.data.frame(gam.fit.covariate(gam.data = gam_df,  
                                                                          region = as.character(x), 
                                                                          smooth_var = smooth.var, 
                                                                          covariate.interest = covariate.interest, 
                                                                          covariates.noninterest = covariates.noninterest, 
                                                                          knots = k, 
                                                                          set_fx = set.fx)) %>% mutate(region = x)}, mc.cores = 4) 
  sex.maineffects <- do.call(rbind, sex.maineffects)
  
  write.csv(sex.maineffects, sprintf("%1$s/%2$s_GAM_%3$s_maineffects_%4$s.csv", GAM_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
  
}

# 2) run gam.factorsmooth.interaction to look at the sex-by-age interaction 
run_gam.factorsmooth.interaction <- function(gam_df, smooth.var, int.var, covs, k, set.fx, filename, atlas){
  agebysex.interactioneffects <- mclapply(parcel.labels, 
                                          function(x){as.data.frame(gam.factorsmooth.interaction(gam.data = gam_df,  
                                                                                                 region = as.character(x), 
                                                                                                 smooth_var = smooth.var, 
                                                                                                 int_var = int.var, 
                                                                                                 covariates = covs, 
                                                                                                 knots = k, 
                                                                                                 set_fx = set.fx)) %>% mutate(region = x)}, mc.cores = 4) 
  agebysex.interactioneffects <- do.call(rbind, agebysex.interactioneffects)
  write.csv(agebysex.interactioneffects, sprintf("%1$s/%2$s_GAM_ageby%3$s_interaction_%4$s.csv", GAM_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
}


# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectory for female and males
run_gam.smooth.predict.sexinteraction <- function(gam_df, smooth_var, int_var, covs, k, set_fx, filename, atlas, age1, age2){
  sex_levels <- c("Male", "Female")   
  
  results_list <- list()
  for (sex in sex_levels) {
    results <- lapply(parcel.labels, function(x) {
      as.data.frame(gam.smooth.predict.covariateinteraction.factor(
        gam.data = gam_df, 
        region = as.character(x), 
        smooth_var = smooth_var, 
        int_var = int_var, 
        int_var.predict = sex,
        covariates = covs, 
        knots = k, 
        set_fx = set_fx,
        increments = 200,
        age1 = age1,
        age2 = age2)) %>% mutate(region = x, sex = sex)
    })
    
    results <- do.call(rbind, results)
    results_list[[sex]] <- results
  }
  
  all_results <- do.call(rbind, results_list)
  write.csv(all_results, sprintf("%1$s/%2$s_GAM_ageby%3$s_fitted_%4$s.csv", GAM_outputs_dir, dataset, filename, atlas), quote = F, row.names =F)
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
demographics <- read.csv(config_data$demographics, stringsAsFactors = TRUE)
demographics <- demographics %>% select(sub, sex, age, meanFD_avgSes)

# load functional connectivity data (accomodating different file names for datasets that got harmonized)
covbat_datasets <- c("HCPD", "HBN")

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
write.table(gam_df$sub, paste0(config_data$covariate_input_root, "/sample_info/final_sex_diff_sample.txt"), col.names=F, row.names=F, quote = TRUE)

# check for sex diff in meanFD_avgSes (average average mean_fd across concatenated sessions)
print(t.test(gam_df$meanFD_avgSes[gam_df$sex=="Male"], gam_df$meanFD_avgSes[gam_df$sex=="Female"]))

smooth.var = "age"
k = 3
set.fx = T
covs = "meanFD_avgSes + sex"

if (dataset == "HBN") {
  age1 = 5
  age2 = 22
} else if (dataset == "HCPD") {
  age1 = 8
  age2 = 22
} else if (dataset == "PNC") {
  age1 = 8
  age2 = 23
} else if (dataset == "NKI") {
  age1 = 6
  age2 = 22
}

################### 
# Fit da GAMz 
###################
# set filename variable for saving out GAM results
filename = "sex"
int.var = "sex" # sex variable already ordered and formatted in ~/covariate_analyses/sex_diff/code/setup/prep_demographic_info.R

# 1) run gam.fit.covariate to look at the main effect of sex  
covariate.interest = "sex"
covariates.noninterest = "meanFD_avgSes"
run_gam.fit.covariate(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx, filename, atlas)

# 2) run gam.factorsmooth.interaction to look at the sex-by-age interaction 
run_gam.factorsmooth.interaction(gam_df, smooth.var, int.var, covs, k, set.fx, filename, atlas)

# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectories for female and males
covs = "sex + meanFD_avgSes"
run_gam.smooth.predict.sexinteraction(gam_df, smooth.var, int.var, covs, k, set.fx, filename, atlas, age1, age2)

print("Script finished!")


