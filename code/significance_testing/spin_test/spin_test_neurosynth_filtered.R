library(mgcv)
library(dplyr)
library(ggplot2)
library(purrr)
library(rjson)
library(stringr)
library(tidyr)
source("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/results/sex_diff_main_figures.R")
source("/cbica/projects/network_replication/software/perm.sphere.p.R")

print("Running spin test for neurosynth")

################## 
# Set Directories 
################## 
proj_root <- "/cbica/projects/network_replication/covariate_analyses/sex_diff/"
input_root <- paste0(proj_root, "input")
output_root <- paste0(proj_root, "output")

output_dir <- "/cbica/projects/network_replication/covariate_analyses/sex_diff/output/all_datasets"
spin_test_outputs_dir <- paste0(output_dir, "/spin_test/neurosynth")

if (dir.exists(spin_test_outputs_dir)) {
  print(paste(spin_test_outputs_dir, "already exists"))
} else {
  dir.create(spin_test_outputs_dir, recursive = TRUE)
  print(paste(spin_test_outputs_dir, "created"))
}

################### 
# Define functions 
###################

sexdiffmap_decoding <- function(df, metric, term, fc_measure){
  
  x <- df[, term]
  y <- df[, metric]
  
  # --- Threshold Neurosynth z-scores ---
  # Keep only z > 1.64 (≈p<0.05 one-tailed)
  x[x < 1.64] <- 0
  
  if (length(x) != nrow(perm.id.full_schaefer200x17)) {
    message(paste("Length mismatch for", term, ":", length(x), "vs rotation", nrow(perm.id.full_schaefer200x17)))
    return(data.frame(term = term, term.correlation = NA, p_spin = NA))
  }
  
  # Safeguard: skip terms with too few suprathreshold parcels
  if (sum(x > 0, na.rm = TRUE) < 5) {
    message(paste("Skipping term:", term, "(too few suprathreshold parcels)"))
    return(data.frame(term = term, term.correlation = NA, p_spin = NA))
  }
  
  # Compute correlation safely
  term.cor <- suppressWarnings(cor(x, y, use = "complete.obs", method = "spearman"))
  
  # Skip spin if correlation invalid or variance zero
  if (is.na(term.cor) || length(unique(x)) < 2 || length(unique(y)) < 2) {
    message(paste("Skipping term:", term, "(no variance or NA correlation)"))
    return(data.frame(term = term, term.correlation = NA, p_spin = NA))
  }
  
  print(paste("Spinning", fc_measure, term))
  
  # Run spin test safely
  term.pspin <- tryCatch(
    perm.sphere.p(x, y, perm.id.full_schaefer200x17, corr.type = "spearman"),
    error = function(e) {
      message(paste("Spin test failed for term:", term, "→", e$message))
      NA
    }
  )
  
  data.frame(
    term = term,
    term.correlation = term.cor,
    p_spin = term.pspin
  )
}

################### 
# Load files
###################

# --- Load and threshold Neurosynth terms ---
print("Loading Neurosynth term z-maps and applying threshold (z < 1.64 set to 0)...")
neurosynth.terms <- read.csv("/cbica/projects/network_replication/covariate_analyses/sex_diff/input/neurosynth_annotations/schaefer200/schaefer200x17_neurosynth_125terms.csv") %>%
  select(-regionID) %>%
  rename(region = label) %>%
  mutate(across(-region, ~ ifelse(. < 1.64, 0, .)))  # threshold here

neurosynth.termlist <- names(neurosynth.terms)[-1]

# Load rotated Schaefer parcellation spins
perm.id.full_schaefer200x17 <- readRDS(
  sprintf("/cbica/projects/network_replication/software/rotate_parcellation/%s.coords_sphericalrotations_N10k_seed10.rds", "schaefer200x17")
)


# load BNC and WNC sex effects 
sex.main.PNC.BNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "PNC", "BNC")) 
sex.main.PNC.BNC$region <- paste0("17", sex.main.PNC.BNC$region) 
sex.main.PNC.BNC <- sig.effects(sex.main.PNC.BNC, "BNC", "sex", "PNC") 

sex.main.HCPD.BNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "HCPD", "BNC")) 
sex.main.HCPD.BNC <- sig.effects(sex.main.HCPD.BNC, "BNC", "sex", "HCPD") 
sex.main.HCPD.BNC$region <- paste0("17", sex.main.HCPD.BNC$region) 

sex.main.NKI.BNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "NKI", "BNC")) 
sex.main.NKI.BNC <- sig.effects(sex.main.NKI.BNC, "BNC", "sex", "NKI") 
sex.main.NKI.BNC$region <- paste0("17", sex.main.NKI.BNC$region) 

sex.main.HBN.BNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "HBN", "BNC")) 
sex.main.HBN.BNC <- sig.effects(sex.main.HBN.BNC, "BNC", "sex", "HBN") 
sex.main.HBN.BNC$region <- paste0("17", sex.main.HBN.BNC$region) 

sex.main.all_datasets.BNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "all_datasets", "BNC")) 
sex.main.all_datasets.BNC <- sig.effects(sex.main.all_datasets.BNC, "BNC", "sex", "all_datasets") 
sex.main.all_datasets.BNC$region <- paste0("17", sex.main.all_datasets.BNC$region) 

sex.main.PNC.WNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "PNC", "WNC")) 
sex.main.PNC.WNC$region <- paste0("17", sex.main.PNC.WNC$region) 
sex.main.PNC.WNC <- sig.effects(sex.main.PNC.WNC, "WNC", "sex", "PNC") 

sex.main.HCPD.WNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "HCPD", "WNC")) 
sex.main.HCPD.WNC <- sig.effects(sex.main.HCPD.WNC, "WNC", "sex", "HCPD") 
sex.main.HCPD.WNC$region <- paste0("17", sex.main.HCPD.WNC$region) 

sex.main.NKI.WNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "NKI", "WNC")) 
sex.main.NKI.WNC <- sig.effects(sex.main.NKI.WNC, "WNC", "sex", "NKI") 
sex.main.NKI.WNC$region <- paste0("17", sex.main.NKI.WNC$region) 

sex.main.HBN.WNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "HBN", "WNC")) 
sex.main.HBN.WNC <- sig.effects(sex.main.HBN.WNC, "WNC", "sex", "HBN") 
sex.main.HBN.WNC$region <- paste0("17", sex.main.HBN.WNC$region) 

sex.main.all_datasets.WNC <- read.csv(sprintf("%1$s/%2$s/GAM/%3$s/%2$s_GAM_sex_maineffects_schaefer200x17.csv", output_root, "all_datasets", "WNC")) 
ex.main.all_datasets.WNC <- sig.effects(sex.main.all_datasets.WNC, "WNC", "sex", "all_datasets") 
sex.main.all_datasets.WNC$region <- paste0("17", sex.main.all_datasets.WNC$region)
 
###################################### 
# Spin and save out - Overlap
###################################### 
BNC_sexeffect_overlap <- identify_sexeffect_overlap("BNC")
WNC_sexeffect_overlap <- identify_sexeffect_overlap("WNC")

neurosynth_BNC <- merge(BNC_sexeffect_overlap, neurosynth.terms, by = "region")
neurosynth_WNC <- merge(WNC_sexeffect_overlap, neurosynth.terms, by = "region")

sexeffect.neurosynth.BNC <- map_dfr(neurosynth.termlist, function(t){
  sexdiffmap_decoding(df = neurosynth_BNC, metric = "mean_t_sig2", term = t, fc_measure = "BNC")}) %>% 
  arrange(desc(term.correlation))

sexeffect.neurosynth.WNC <- map_dfr(neurosynth.termlist, function(t){
  sexdiffmap_decoding(df = neurosynth_WNC, metric = "mean_t_sig2", term = t, fc_measure = "WNC")}) %>% 
  arrange(desc(term.correlation))

write.csv(sexeffect.neurosynth.BNC, sprintf("%s/BNC_neurosynth_spin_test_results_overlap_filtered.csv", spin_test_outputs_dir))
write.csv(sexeffect.neurosynth.WNC, sprintf("%s/WNC_neurosynth_spin_test_results_overlap_filtered.csv", spin_test_outputs_dir))

###################################### 
# Spin and save out - Pooled
###################################### 
neurosynth_BNC_pooled <- merge(sex.main.all_datasets.BNC, neurosynth.terms, by = "region")
neurosynth_WNC_pooled <- merge(sex.main.all_datasets.WNC, neurosynth.terms, by = "region")

sexeffect.neurosynth.BNC.pooled <- map_dfr(neurosynth.termlist, function(t){
  sexdiffmap_decoding(df = neurosynth_BNC_pooled, metric = "GAM.cov.tvalue", term = t, fc_measure = "BNC")}) %>% 
  arrange(desc(term.correlation))

sexeffect.neurosynth.WNC.pooled <- map_dfr(neurosynth.termlist, function(t){
  sexdiffmap_decoding(df = neurosynth_WNC_pooled, metric = "GAM.cov.tvalue", term = t, fc_measure = "WNC")}) %>% 
  arrange(desc(term.correlation))

write.csv(sexeffect.neurosynth.BNC.pooled, sprintf("%s/BNC_neurosynth_spin_test_results_pooled_filtered.csv", spin_test_outputs_dir))
write.csv(sexeffect.neurosynth.WNC.pooled, sprintf("%s/WNC_neurosynth_spin_test_results_pooled_filtered.csv", spin_test_outputs_dir))
