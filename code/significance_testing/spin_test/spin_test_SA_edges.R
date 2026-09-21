library(mgcv)
library(dplyr)
library(ggplot2)
library(rjson)
library(stringr)
library(tidyr)
source("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/results/sex_diff_main_figures.R")

# a script to run spin test How much of the variance in 
# edgewise sex effects is captured by that smooth hierarchical organization -- 
# or how strongly the S–A axis explains the spatial distribution of sex differences in edge-level connectivity

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
atlas = args[2]
metric = "edge" 

print(paste("Running spin test for", dataset))

################## 
# Set Directories 
################## 
if (dataset != "all_datasets") {
  config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
  outputs_root <- config_data$covariate_output_root
  
} else if (dataset == "all_datasets") {
  outputs_root <- "/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets"
}

spin_test_outputs_dir <- paste0(outputs_root, "/spin_test/", metric)
metric_outputs_dir <- paste0(outputs_root, "/", metric)

if (dir.exists(spin_test_outputs_dir)) {
  print(paste(spin_test_outputs_dir, "already exists"))
} else {
  dir.create(spin_test_outputs_dir, recursive = TRUE)
  print(paste(spin_test_outputs_dir, "created"))
}

################### 
# Define functions 
###################
# edge level analysis: 
# make a double df for symmetry - parcel1 and parcel2.SA.vec repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
make_double_df <- function(df, fc_measure, covariate, dataset) {
  age_SA.diff <- cbind(df, SA.diff)
  age_SA.diff <- sig.effects(age_SA.diff, fc_measure, covariate, dataset) 
  
  # make a double df for symmetry - parcel1 and parcel2.SA.vec repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
  age_SA.diff2 <- age_SA.diff 
  age_SA.diff2$parcel1.SA.vec <- age_SA.diff$parcel2.SA.vec
  age_SA.diff2$parcel2.SA.vec <- age_SA.diff$parcel1.SA.vec
  
  # "stacked" df
  double_age_SA.diff <- rbind(age_SA.diff2,age_SA.diff)
  
  return(double_age_SA.diff)
}


# spin test
run_spin_test <- function(double_df, SA.vec, perm.id.full,
                          measure = "GAM.cov.tvalue", k_basis = 3) {
  p1_idx <- match(double_df$parcel1_id, seq_along(SA.vec))
  p2_idx <- match(double_df$parcel2_id, seq_along(SA.vec))
  p1_SA  <- SA.vec[p1_idx]
  p2_SA  <- SA.vec[p2_idx]
  y <- double_df[[measure]]
  
  g2_obs <- gam(y ~ te(p2_SA, p1_SA, k = k_basis), method = "REML")
  R2_obs <- summary(g2_obs)$r.sq
  
  K <- ncol(perm.id.full)
  R2_null <- numeric(K)
  pb <- txtProgressBar(min = 0, max = K, style = 3)
  message("Running spin test (R²) with ", K, " spins ...")
  
  for (k in seq_len(K)) {
    perm_idx <- perm.id.full[, k]
    SA_spin <- SA.vec[perm_idx]
    p1_SA_k <- SA_spin[p1_idx]
    p2_SA_k <- SA_spin[p2_idx]
    g2_k <- try(gam(y ~ te(p2_SA_k, p1_SA_k, k = k_basis), method = "REML"), silent = TRUE)
    if (!inherits(g2_k, "try-error")) {
      R2_null[k] <- summary(g2_k)$r.sq
    } else {
      R2_null[k] <- NA
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  valid_R2 <- is.finite(R2_null)
  p_R2 <- (1 + sum(R2_null[valid_R2] >= R2_obs)) / (sum(valid_R2) + 1)
  
  cat("\nObserved stats:\n",
      sprintf("  R² = %.4f (p_spin = %.4f)\n", R2_obs, p_R2),
      sprintf("\nConverged spins: %d / %d\n", sum(valid_R2), K))
  
  return(list(
    observed = list(R2 = R2_obs),
    p_values = list(R2 = p_R2),
    nulls = list(R2 = R2_null)
  ))
}




################### 
# Load files
###################

# load S-A axis
SAaxis <- read.csv(sprintf("/cbica/projects/network_replication/SAaxis/%s_SAaxis.csv", atlas))
SAaxis <- SAaxis %>% rename(region=label) %>% rename(parcel_id = X)

# load rotated schaefer200x17
perm.id.full <- readRDS(sprintf("/cbica/projects/network_replication/software/rotate_parcellation/%s.coords_sphericalrotations_N10k_seed10.rds", atlas))

# load delta-SA rank for each edge
SA.diff <- read.csv(sprintf("/cbica/projects/network_replication/atlases/edge/SA.diff_%s.csv", atlas))

# load edges
if(dataset != "all_datasets") {
  sex.main.edge <- read.csv(sprintf("%1$s/GAM/%3$s/%2$s_GAM_sex_maineffects.csv", outputs_root, dataset, metric))
} else if (dataset == "all_datasets") {
  sex.main.edge <- read.csv(sprintf("%1$s/GAM/%3$s/%2$s_GAM_sex_maineffects_%4$s.csv", outputs_root, dataset, metric, atlas))
}

# make double df for symmetry 
double_df <- make_double_df(df = sex.main.edge, fc_measure = metric, covariate = "sex", dataset = dataset)
 
################### 
# Spin and save out
###################
spin <- run_spin_test(double_df = double_df, SA.vec = SAaxis$SA.axis_rank, perm.id.full)
 
saveRDS(spin, sprintf("%s/SA_edges_%s_spin_test_results.RData", spin_test_outputs_dir, atlas))
