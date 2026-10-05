library(ComBatFamily)
library(data.table)
library(dplyr)
library(mgcv)
library(rjson)
library(stringr)
library(tidyr)

##################
# Set variables
##################
# Usage:
#   Rscript covbat_fc_matrices.R HCPD schaefer200x17
#   Rscript covbat_fc_matrices.R HBN  schaefer200x17
#   Rscript covbat_fc_matrices.R all_datasets schaefer200x17
#
# NOTE: PNC and NKI are single-site and do not need this script --
# only run for datasets/pooling where a site variable exists.

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]          # "HCPD", "HBN", or "all_datasets"
atlas <- args[2]            # e.g. "schaefer200x17"
n_parcels <- as.integer(str_extract(atlas, "(?<=schaefer)[0-9]+"))
print(paste("Harmonizing FC matrices for", dataset, atlas))

##################
# Helpers
##################

# Flatten a subject's upper triangle (no diagonal) into a named vector.
# Matches the ordering used by extractParcel2ParcelConn() elsewhere in
# this project, so edge names/order are consistent across analyses.
matrix_to_edge_vector <- function(m, edge_names) {
  v <- m[upper.tri(m)]
  names(v) <- edge_names
  v
}

# Inverse: rebuild a symmetric matrix (zero diagonal) from an edge vector.
edge_vector_to_matrix <- function(v, n) {
  m <- matrix(0, n, n)
  m[upper.tri(m)] <- v
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

load_matrices_for <- function(config, key) {
  # Returns list(subxedge_uppertri = data.frame with 'subject' + edge cols,
  #               ids = subject id vector, n = parcel count, edge_names)
  #
  # NOTE: "subxedge_uppertri" is deliberately NOT called "subxedge" --
  # that name is already used elsewhere in this project for the
  # lower-triangle, value-matched edge table produced by
  # extractParcel2ParcelConn() (see compute_connectivity_metrics.R). This
  # object uses upper.tri() with explicit index-based exclusion instead,
  # so the two are NOT column-for-column interchangeable. Keeping the
  # names distinct avoids a silent overwrite/mismatch if both pipelines
  # are ever sourced in the same R session.
  conn_dir <- file.path(config$conn_matrices_dir)
  sample_ids <- readLines(config$subject_list)
  sample_ids <- gsub('"', '', trimws(sample_ids))
  sample_ids <- sample_ids[grepl("[0-9]", sample_ids)]
  if (config$dataset == "NKI") { # I annoyingly named the NKI connectivity matrices in a different format back in the day lol
    sample_ids <- sub("^sub-", "", sample_ids)
  }

  edge_names <- NULL
  rows <- list()
  used <- character(0)

  
  for (sid in sample_ids) {
    p <- sprintf("%s/%s_ConnMatrices.RData", conn_dir, sid)
    if (!file.exists(p)) next
    
    if (config$dataset == "NKI") { # naming convention inside NKI's matrices is also different from the other datasets due to presence of sessions
      x <- readRDS(p)
      match <- grep(key, names(x), value = TRUE, fixed = TRUE)
      if (length(match) != 1) next
      m <- x[[match]]
    } else {
      m <- readRDS(p)[[key]]
    }
    
    if (is.null(edge_names)) {
      combos <- which(upper.tri(m), arr.ind = TRUE)
      edge_names <- paste0("e", combos[, 1], "_", combos[, 2])
    }
    
    rows[[length(rows) + 1]] <- matrix_to_edge_vector(m, edge_names)
    used <- c(used, sid)
  }
 

  subxedge_uppertri <- as.data.frame(do.call(rbind, rows))
  subxedge_uppertri$subject <- used
  subxedge_uppertri <- subxedge_uppertri %>% relocate(subject)
  list(subxedge_uppertri = subxedge_uppertri, ids = used, n = nrow(m), edge_names = edge_names)
}

# key in the RData list, matching your existing 17-network convention
# (schaefer200 -> schaefer217_conn, schaefer400 -> schaefer417_conn)
atlas_key <- function(atlas) {
  n <- str_extract(atlas, "schaefer[0-9]")
  paste0(n, "17_conn")
}

run_covbat <- function(subxedge_uppertri, demographics, batch_vec) {
  subxedge_uppertri <- subxedge_uppertri[complete.cases(subxedge_uppertri), ]

  if (!identical(as.character(demographics$sub), as.character(subxedge_uppertri$subject))) {
    stop("ERROR: subject IDs in demographics and subxedge_uppertri do not match. Exiting.")
  }

  edge_mat <- data.frame(subxedge_uppertri)[, -1]
  row.names(edge_mat) <- demographics$sub

  covar_df <- data.frame(
    sub = demographics$sub,
    age = as.numeric(demographics$age),
    sex = as.factor(demographics$sex),
    meanFD_avgSes = as.numeric(demographics$meanFD_avgSes)
  )

  print(sprintf("running covfam on %d subjects x %d edges", nrow(edge_mat), ncol(edge_mat)))
  data.harmonized <- covfam(edge_mat, bat = as.factor(batch_vec), covar = covar_df,
                            model = gam,
                            formula = y ~ s(age, k = 3, fx = TRUE) + as.factor(sex) + as.numeric(meanFD_avgSes))
  print("covbat complete")

  out <- data.frame(data.harmonized$dat.covbat)
  out$subject <- rownames(edge_mat)
  out <- out %>% relocate(subject)
  rownames(out) <- NULL
  out
}


##################
# Single-dataset path (HCPD, HBN)
##################
if (dataset != "all_datasets") {

  config_data <- fromJSON(file = sprintf(
    "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_%s.json", dataset))
  outputs_root <- config_data$covariate_output_root
  fc_outputs_dir <- paste0(outputs_root, "/fc_matrices")
  if (!dir.exists(fc_outputs_dir)) dir.create(fc_outputs_dir, recursive = TRUE)

  demographics <- read.csv(config_data$demographics, stringsAsFactors = TRUE)
  if (dataset == "HBN") demographics <- demographics %>% rename(site = ses)
  if (!"site" %in% names(demographics)) {
    stop(sprintf("ERROR: %s has no site column -- this dataset may not need harmonization.", dataset))
  }

  loaded <- load_matrices_for(config_data, atlas_key(atlas))
  n <- loaded$n
  edge_names <- loaded$edge_names

  demographics <- demographics %>% filter(sub %in% loaded$ids) %>%
    arrange(match(sub, loaded$ids))

  covbat_out <- run_covbat(loaded$subxedge_uppertri, demographics, droplevels(demographics$site))

  saveRDS(covbat_out, sprintf("%s/subxedge_uppertri_%s_covbat.RData", fc_outputs_dir, atlas))
  print(paste("saved edge-level covbat output ->", fc_outputs_dir))

  # ---- rebuild per-subject matrices for the gradient pipeline ----
  mats_dir <- file.path(fc_outputs_dir, "matrices_covbat")
  if (!dir.exists(mats_dir)) dir.create(mats_dir, recursive = TRUE)

  for (i in seq_len(nrow(covbat_out))) {
    sid <- covbat_out$subject[i]
    v <- as.numeric(covbat_out[i, edge_names])
    m <- edge_vector_to_matrix(v, n)
    saveRDS(m, sprintf("%s/%s_fc_covbat.RData", mats_dir, sid))
  }
  print(sprintf("saved %d harmonized matrices -> %s", nrow(covbat_out), mats_dir))
}


##################
# Pooled path (all_datasets)
##################
if (dataset == "all_datasets") {

  configs <- list(
    PNC  = fromJSON(file = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_PNC.json"),
    HCPD = fromJSON(file = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_HCPD.json"),
    NKI  = fromJSON(file = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_NKI.json"),
    HBN  = fromJSON(file = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/code/config/config_HBN.json")
  )

  outputs_root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff/output/all_datasets"
  fc_outputs_dir <- paste0(outputs_root, "/fc_matrices")
  if (!dir.exists(fc_outputs_dir)) dir.create(fc_outputs_dir, recursive = TRUE)

  demo_list <- list()
  edge_list <- list()
  n_ref <- NULL
  edge_names_ref <- NULL

  for (ds in names(configs)) {
    cfg <- configs[[ds]]
    demo <- read.csv(cfg$demographics, stringsAsFactors = TRUE)
    if (ds == "HBN") demo <- demo %>% rename(site = ses)
    # single-site datasets get their own dataset name as the batch label,
    # matching the dataset_site convention in covbat_all_datasets.R
    demo$dataset_site <- if ("site" %in% names(demo)) paste0(ds, "_", demo$site) else ds
    demo <- demo %>% select(sub, age, sex, meanFD_avgSes, dataset_site)

    loaded <- load_matrices_for(cfg, atlas_key(atlas))
    if (is.null(n_ref)) { n_ref <- loaded$n; edge_names_ref <- loaded$edge_names }
    stopifnot(loaded$n == n_ref)   # all datasets must share the same parcellation

    demo <- demo %>% filter(sub %in% loaded$ids) %>% arrange(match(sub, loaded$ids))

    demo_list[[ds]] <- demo
    edge_list[[ds]] <- loaded$subxedge_uppertri
    print(sprintf("%s: %d subjects loaded", ds, nrow(demo)))
  }

  all_demographics <- bind_rows(demo_list)
  all_subxedge_uppertri <- bind_rows(edge_list)
  
  cat("dataset_site counts going into covbat:\n")
  print(table(all_demographics$dataset_site))

  covbat_out <- run_covbat(all_subxedge_uppertri, all_demographics, droplevels(as.factor(all_demographics$dataset_site)))
  saveRDS(covbat_out, sprintf("%s/subxedge_uppertri_%s_covbat.RData", fc_outputs_dir, atlas))
  print(paste("saved pooled edge-level covbat output ->", fc_outputs_dir))

  mats_dir <- file.path(fc_outputs_dir, "matrices_covbat")
  if (!dir.exists(mats_dir)) dir.create(mats_dir, recursive = TRUE)

  for (i in seq_len(nrow(covbat_out))) {
    sid <- covbat_out$subject[i]
    v <- as.numeric(covbat_out[i, edge_names_ref])
    m <- edge_vector_to_matrix(v, n_ref)
    saveRDS(m, sprintf("%s/%s_fc_covbat.RData", mats_dir, sid))
  }
  print(sprintf("saved %d harmonized matrices -> %s", nrow(covbat_out), mats_dir))
}
