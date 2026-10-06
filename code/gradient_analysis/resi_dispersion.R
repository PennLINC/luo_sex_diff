library(dplyr)
library(tidyr)
library(mgcv)
library(RESI)
library(parallel)
source("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/compute_effect_size/resiPEse_Th1_functions.R")

##################
# This script computes RESI effect size for global, within- and between-network dispersion
# OUTPUT:
# - <DATASET>_resi_sex_dispersion_global.csv: RESI sex effect and associated statistics for global dispersion.
# - <DATASET>_resi_dispersion_meanfc_comparison.csv: Comparison of the global-dispersion sex effect with versus without adjustment for mean FC.
# - <DATASET>_resi_sex_dispersion_within_network.csv: RESI sex effects for each within-network dispersion measure, including Bonferroni and FDR-adjusted p-values.
# - <DATASET>_resi_sex_dispersion_between_network.csv: RESI sex effects for each pairwise between-network dispersion measure, including Bonferroni and FDR-adjusted p-values.
##################


##################
# Set variables
##################
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]        # PNC, NKI, HCPD, HBN, all_datasets

outputs_root <- sprintf("/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion", dataset)
resi_outputs_dir <- paste0(outputs_root, "/resi")
if (!dir.exists(resi_outputs_dir)) dir.create(resi_outputs_dir, recursive = TRUE)

k <- 3
set.fx <- TRUE
covariate.interest <- "sex"

# PRIMARY covariate set excludes mean_fc
covariates.primary <- "meanFD_avgSes"

# sensitivity model additionally adjusts for mean_fc
covariates.with_meanfc <- "meanFD_avgSes + mean_fc"

##################
# Shared RESI helper
##################
# Positive RESI = higher in females; negative RESI = higher in males.
run_resi <- function(outcome_var, data, covariates.noninterest = covariates.primary) {
  formula <- as.formula(paste0("`", outcome_var, "` ~ s(age, k = ", k, ", fx = ", set.fx,
                               ") + ", covariates.noninterest, " + ", covariate.interest))
  mod <- mgcv::gam(formula, data = data, method = "REML")
  
  th1 <- resiPEse_Th1(mod, variable = "sexFemale", unsigned = FALSE, torz = "t", type = "HC0")
  
  t_robust <- th1$Estimate / th1$Std.Error
  p_robust <- 2 * pt(-abs(t_robust), df = mod$df.residual)
  t2s <- RESI::t2S(t_robust, rdf = mod$df.residual, n = nrow(data), unbiased = TRUE)
  
  data.frame(outcome = outcome_var, covariates = covariates.noninterest,
             RESI_Xinyu = th1$RESI_Xinyu, LCI = th1$LCI, UCI = th1$UCI,
             RESI_t2S = as.numeric(t2s), beta = th1$Estimate, se = th1$Std.Error,
             t = t_robust, p = p_robust)
}

##################
# PART 1 -- Global dispersion
##################
gam_df <- read.csv(sprintf("%s/dispersion_global_%s.csv", outputs_root, dataset))
gam_df <- gam_df %>% drop_na(sex)
gam_df$sex <- factor(gam_df$sex, levels = c("Male", "Female"))

cat(sprintf("\n=== %s: motion sanity check ===\n", dataset))
print(t.test(gam_df$meanFD_avgSes[gam_df$sex == "Male"],
             gam_df$meanFD_avgSes[gam_df$sex == "Female"]))

cat("\n=== GLOBAL DISPERSION: PRIMARY MODEL ===\n")
global_result <- run_resi("dispersion", gam_df)
print(global_result)

write.csv(global_result,
          sprintf("%s/%s_resi_sex_dispersion_global.csv", resi_outputs_dir, dataset),
          quote = FALSE, row.names = FALSE)

##################
# PART 1b -- Global dispersion sensitivity: with vs without mean_fc
##################
cat("\n=== GLOBAL DISPERSION: with vs without mean_fc ===\n")

res_no_meanfc <- run_resi("dispersion", gam_df, covariates.noninterest = covariates.primary)
res_with_meanfc <- run_resi("dispersion", gam_df, covariates.noninterest = covariates.with_meanfc)

meanfc_compare <- rbind(
  data.frame(model = "without mean_fc", res_no_meanfc[c("RESI_Xinyu", "beta", "p")]),
  data.frame(model = "with mean_fc", res_with_meanfc[c("RESI_Xinyu", "beta", "p")])
)

print(meanfc_compare, row.names = FALSE)

pct_change <- (abs(res_with_meanfc$beta) - abs(res_no_meanfc$beta)) / abs(res_no_meanfc$beta) * 100
cat(sprintf("\n|beta| changed %+.1f%% when adding mean_fc (%s)\n", pct_change,
            if (pct_change < -10) "ATTENUATED -- possible confound"
            else if (pct_change > 10) "STRENGTHENED"
            else "STABLE"))

write.csv(meanfc_compare,
          sprintf("%s/%s_resi_dispersion_meanfc_comparison.csv", resi_outputs_dir, dataset),
          quote = FALSE, row.names = FALSE)

##################
# PART 2 -- Within- and between-network dispersion
# PRIMARY models exclude mean_fc
##################
net_df <- read.csv(sprintf("%s/dispersion_network_%s.csv", outputs_root, dataset), check.names = FALSE)
net_df <- net_df %>% drop_na(sex)
net_df$sex <- factor(net_df$sex, levels = c("Male", "Female"))

wn_cols <- grep("^wn_", names(net_df), value = TRUE)
bn_cols <- grep("^bn_", names(net_df), value = TRUE)

cat(sprintf("\n%d within-network outcomes, %d between-network pairs, n = %d\n",
            length(wn_cols), length(bn_cols), nrow(net_df)))

n_cores <- max(1, detectCores() - 2)

cat("\n=== WITHIN-NETWORK DISPERSION: PRIMARY MODEL ===\n")
wn_results <- do.call(rbind, mclapply(wn_cols, run_resi, data = net_df, mc.cores = n_cores))
wn_results$p_bonf <- p.adjust(wn_results$p, method = "bonferroni")
wn_results$p_fdr <- p.adjust(wn_results$p, method = "fdr")
wn_results <- wn_results[order(wn_results$p), ]
print(wn_results, row.names = FALSE)

cat("\n=== BETWEEN-NETWORK DISPERSION: PRIMARY MODEL ===\n")
bn_results <- do.call(rbind, mclapply(bn_cols, run_resi, data = net_df, mc.cores = n_cores))
bn_results$p_bonf <- p.adjust(bn_results$p, method = "bonferroni")
bn_results$p_fdr <- p.adjust(bn_results$p, method = "fdr")
bn_results <- bn_results[order(bn_results$p), ]
print(head(bn_results, 15), row.names = FALSE)

write.csv(wn_results,
          sprintf("%s/%s_resi_sex_dispersion_within_network.csv", resi_outputs_dir, dataset),
          quote = FALSE, row.names = FALSE)

write.csv(bn_results,
          sprintf("%s/%s_resi_sex_dispersion_between_network.csv", resi_outputs_dir, dataset),
          quote = FALSE, row.names = FALSE)

cat(sprintf("\n%s done. global RESI (raw dispersion, no mean_fc) = %+.4f (p=%.4g)\n",
            dataset, global_result$RESI_Xinyu, global_result$p))

cat(sprintf("within-network: %d/%d survive FDR\n",
            sum(wn_results$p_fdr < 0.05), nrow(wn_results)))

cat(sprintf("between-network: %d/%d survive FDR\n",
            sum(bn_results$p_fdr < 0.05), nrow(bn_results)))