library(mgcv)
library(dplyr)
library(purrr)
library(rjson)
library(stringr)
library(tidyr)

##################
# Set directories
##################

# Set project root depending on whether code is running on PARCC or locally
if (dir.exists("/ceph/projects/sattertt/pennlinc-parcc/network_replication")) {
  root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication"
} else if (dir.exists("~/parcc/network_replication")) {
  root <- path.expand("~/parcc/network_replication")
} else {
  stop("Could not find network_replication project.")
}

proj_root <- file.path(root, "covariate_analyses", "sex_diff")
input_root <- file.path(proj_root, "input")
output_root <- file.path(proj_root, "output")
output_dir <- file.path(output_root, "all_datasets")
spin_test_outputs_dir <- file.path(output_dir, "spin_test", "neurosynth")

source(file.path(proj_root, "code", "results", "main_figures.R"))
source(file.path(root, "software", "perm.sphere.p.R"))

if (!dir.exists(spin_test_outputs_dir)) dir.create(spin_test_outputs_dir, recursive = TRUE)

print("Running spin test for neurosynth")

###################
# Define functions
###################
  
# function to correlate each Neurosynth term map with the regional sex-effect map and use
# spin testing to test for significance
# note: FDR correction happens later!
  sexdiffmap_decoding <- function(df, metric, term, fc_measure) { 
    term.cor <- cor(df[, term], df[, metric], use = "complete.obs", method = "spearman")
    print(paste("Spinning", fc_measure, term)) 
    term.pspin <- perm.sphere.p(df[, term], df[, metric], perm.id.full_schaefer200x17, corr.type = "spearman")
    term.results <- data.frame(term = term, term.correlation = term.cor, p_spin = term.pspin)
    return(term.results) 
  }

###################
# Load files
###################

print("Loading Neurosynth term z-maps")
neurosynth.terms <- read.csv(file.path(input_root, "neurosynth_annotations", "schaefer200", "schaefer200x17_neurosynth_125terms.csv")) %>%
  select(-regionID) %>%
  rename(region = label)

neurosynth.termlist <- names(neurosynth.terms)[-1]

# Load rotated Schaefer parcellation spins
perm.id.full_schaefer200x17 <- readRDS(file.path(proj_root, "software", "rotate_parcellation", "schaefer200x17.coords_sphericalrotations_N10k_seed10.rds"))

# Load BNC sex effects
sex.main.PNC.BNC <- read.csv(file.path(output_root, "PNC", "GAM", "BNC", "PNC_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.PNC.BNC$region <- paste0("17", sex.main.PNC.BNC$region)
sex.main.PNC.BNC <- sig.effects(sex.main.PNC.BNC, "BNC", "sex", "PNC")
resi.PNC.BNC <- read.csv(file.path(output_root, "PNC", "resi", "BNC", "PNC_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.PNC.BNC$region <- paste0("17", resi.PNC.BNC$region)
resi.PNC.BNC <- resi.PNC.BNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.PNC.BNC <- merge(sex.main.PNC.BNC, resi.PNC.BNC)

sex.main.HCPD.BNC <- read.csv(file.path(output_root, "HCPD", "GAM", "BNC", "HCPD_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.HCPD.BNC <- sig.effects(sex.main.HCPD.BNC, "BNC", "sex", "HCPD")
sex.main.HCPD.BNC$region <- paste0("17", sex.main.HCPD.BNC$region)
resi.HCPD.BNC <- read.csv(file.path(output_root, "HCPD", "resi", "BNC", "HCPD_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.HCPD.BNC$region <- paste0("17", resi.HCPD.BNC$region)
resi.HCPD.BNC <- resi.HCPD.BNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.HCPD.BNC <- merge(sex.main.HCPD.BNC, resi.HCPD.BNC)

sex.main.NKI.BNC <- read.csv(file.path(output_root, "NKI", "GAM", "BNC", "NKI_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.NKI.BNC <- sig.effects(sex.main.NKI.BNC, "BNC", "sex", "NKI")
sex.main.NKI.BNC$region <- paste0("17", sex.main.NKI.BNC$region)
resi.NKI.BNC <- read.csv(file.path(output_root, "NKI", "resi", "BNC", "NKI_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.NKI.BNC$region <- paste0("17", resi.NKI.BNC$region)
resi.NKI.BNC <- resi.NKI.BNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.NKI.BNC <- merge(sex.main.NKI.BNC, resi.NKI.BNC)

sex.main.HBN.BNC <- read.csv(file.path(output_root, "HBN", "GAM", "BNC", "HBN_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.HBN.BNC <- sig.effects(sex.main.HBN.BNC, "BNC", "sex", "HBN")
sex.main.HBN.BNC$region <- paste0("17", sex.main.HBN.BNC$region)
resi.HBN.BNC <- read.csv(file.path(output_root, "HBN", "resi", "BNC", "HBN_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.HBN.BNC$region <- paste0("17", resi.HBN.BNC$region)
resi.HBN.BNC <- resi.HBN.BNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.HBN.BNC <- merge(sex.main.HBN.BNC, resi.HBN.BNC)
 
# Load WNC sex effects
sex.main.PNC.WNC <- read.csv(file.path(output_root, "PNC", "GAM", "WNC", "PNC_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.PNC.WNC$region <- paste0("17", sex.main.PNC.WNC$region)
sex.main.PNC.WNC <- sig.effects(sex.main.PNC.WNC, "WNC", "sex", "PNC")
resi.PNC.WNC <- read.csv(file.path(output_root, "PNC", "resi", "WNC", "PNC_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.PNC.WNC$region <- paste0("17", resi.PNC.WNC$region)
resi.PNC.WNC <- resi.PNC.WNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.PNC.WNC <- merge(sex.main.PNC.WNC, resi.PNC.WNC)

sex.main.HCPD.WNC <- read.csv(file.path(output_root, "HCPD", "GAM", "WNC", "HCPD_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.HCPD.WNC <- sig.effects(sex.main.HCPD.WNC, "WNC", "sex", "HCPD")
sex.main.HCPD.WNC$region <- paste0("17", sex.main.HCPD.WNC$region)
resi.HCPD.WNC <- read.csv(file.path(output_root, "HCPD", "resi", "WNC", "HCPD_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.HCPD.WNC$region <- paste0("17", resi.HCPD.WNC$region)
resi.HCPD.WNC <- resi.HCPD.WNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.HCPD.WNC <- merge(sex.main.HCPD.WNC, resi.HCPD.WNC)

sex.main.NKI.WNC <- read.csv(file.path(output_root, "NKI", "GAM", "WNC", "NKI_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.NKI.WNC <- sig.effects(sex.main.NKI.WNC, "WNC", "sex", "NKI")
sex.main.NKI.WNC$region <- paste0("17", sex.main.NKI.WNC$region)
resi.NKI.WNC <- read.csv(file.path(output_root, "NKI", "resi", "WNC", "NKI_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.NKI.WNC$region <- paste0("17", resi.NKI.WNC$region)
resi.NKI.WNC <- resi.NKI.WNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.NKI.WNC <- merge(sex.main.NKI.WNC, resi.NKI.WNC)

sex.main.HBN.WNC <- read.csv(file.path(output_root, "HBN", "GAM", "WNC", "HBN_GAM_sex_maineffects_schaefer200x17.csv"))
sex.main.HBN.WNC <- sig.effects(sex.main.HBN.WNC, "WNC", "sex", "HBN")
sex.main.HBN.WNC$region <- paste0("17", sex.main.HBN.WNC$region)
resi.HBN.WNC <- read.csv(file.path(output_root, "HBN", "resi", "WNC", "HBN_resiPEse_Th1_sex_maineffects_schaefer200x17.csv"))
resi.HBN.WNC$region <- paste0("17", resi.HBN.WNC$region)
resi.HBN.WNC <- resi.HBN.WNC %>% rename(RESI = RESI_Xinyu) %>% select(region, RESI)
sex.main.HBN.WNC <- merge(sex.main.HBN.WNC, resi.HBN.WNC)
 
######################################
# Spin and save out - Overlap
######################################

BNC_sexeffect_overlap <- identify_sexeffect_overlap("BNC")
WNC_sexeffect_overlap <- identify_sexeffect_overlap("WNC")

neurosynth_BNC <- merge(BNC_sexeffect_overlap, neurosynth.terms, by = "region")
neurosynth_WNC <- merge(WNC_sexeffect_overlap, neurosynth.terms, by = "region")

sexeffect.neurosynth.BNC <- map_dfr(neurosynth.termlist, function(t) {
  sexdiffmap_decoding(df = neurosynth_BNC, metric = "mean_resi_sig2", term = t, fc_measure = "BNC")
}) %>% arrange(desc(term.correlation))

sexeffect.neurosynth.WNC <- map_dfr(neurosynth.termlist, function(t) {
  sexdiffmap_decoding(df = neurosynth_WNC, metric = "mean_resi_sig2", term = t, fc_measure = "WNC")
}) %>% arrange(desc(term.correlation))

write.csv(sexeffect.neurosynth.BNC, file.path(spin_test_outputs_dir, "BNC_neurosynth_spin_test_results_overlap.csv"), row.names = FALSE)
write.csv(sexeffect.neurosynth.WNC, file.path(spin_test_outputs_dir, "WNC_neurosynth_spin_test_results_overlap.csv"), row.names = FALSE)
 