library(mgcv)
library(dplyr)
library(ggplot2)
library(rjson)
library(stringr)
library(tidyr)
source("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/results/sex_diff_draft_figures.R")

# a script to run spin test how much of the variance in 
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
  config_data <- fromJSON(file=sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
  outputs_root <- config_data$covariate_output_root
  
} else if (dataset == "all_datasets") {
  outputs_root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/output/all_datasets"
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

 
# Spin test of the relationship between S-A axis and edgewise sex effects:
# Fits a GAM relating edgewise sex effects to the S-A rank of both edge endpoints (using te())
# Parcel-level S-A values are spatially rotated and reassigned to the fixed edge
# endpoints to generate a spatially constrained null distribution of explained variance (R-sq).
run_spin_test <- function(double_df, SA.vec, perm.id.full, measure = "GAM.cov.tvalue", k_basis = 3) {
  
  # Match edge endpoints to parcel-level S-A values (assumes SA.vec is ordered by parcel ID)
  p1_idx <- match(double_df$parcel1_id, seq_along(SA.vec))
  p2_idx <- match(double_df$parcel2_id, seq_along(SA.vec))
  stopifnot(!anyNA(p1_idx), !anyNA(p2_idx))
  
  p1_SA <- SA.vec[p1_idx]
  p2_SA <- SA.vec[p2_idx]
  y <- double_df[[measure]]
  
  # Observed model: variance in edgewise sex effects explained by the S-A positions
  # of the two connected parcels
  g2_obs <- gam(y ~ te(p2_SA, p1_SA, k = k_basis), method = "REML")
  R2_obs <- summary(g2_obs)$r.sq
  
  # Generate null distribution by spatially rotating the parcel-level S-A map
  K <- ncol(perm.id.full)
  R2_null <- rep(NA_real_, K)
  pb <- txtProgressBar(min = 0, max = K, style = 3)
  message("Running spin test (R-sq) with ", K, " spins ...")
  
  for (k in seq_len(K)) {
    SA_spin <- SA.vec[perm.id.full[, k]]
    
    # Reassign spun S-A values to the endpoints of each fixed edge
    p1_SA_k <- SA_spin[p1_idx]
    p2_SA_k <- SA_spin[p2_idx]
    
    # Refit the same GAM and retain explained variance
    g2_k <- try(gam(y ~ te(p2_SA_k, p1_SA_k, k = k_basis), method = "REML"), silent = TRUE)
    if (!inherits(g2_k, "try-error")) R2_null[k] <- summary(g2_k)$r.sq
    
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  # Compare observed R² with the spatially constrained null distribution
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
SAaxis <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/SAaxis/%s_SAaxis.csv", atlas))
SAaxis <- SAaxis %>% rename(region=label) %>% rename(parcel_id = X)

# load rotated schaefer200x17
perm.id.full <- readRDS(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/software/rotate_parcellation/%s.coords_sphericalrotations_N10k_seed10.rds", atlas))

# load delta-SA rank for each edge
SA.diff <- read.csv(sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/edge/SA.diff_%s.csv", atlas))

# load edges
if(dataset != "all_datasets") {
  sex.main.edge <- read.csv(sprintf("%1$s/GAM/%3$s/%2$s_GAM_sex_maineffects_%4$s.csv", outputs_root, dataset, metric, atlas))
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
