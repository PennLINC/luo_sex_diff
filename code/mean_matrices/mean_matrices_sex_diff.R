library(data.table)
library(dplyr)
library(rjson)
library(tidyr)

# This script computes mean matrices for each sex.

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
atlas = args[2]
print(paste("Computing sex-specific mean matrices for", dataset, "in", atlas))

################## 
# Set directories  
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/code/config/config_%1$s.json", dataset))
mat_dir <- file.path("/cbica/projects/network_replication/manuscript/input", dataset, "connMatricesData", "connectivity_matrices")
out_dir <- paste0(config_data$covariate_input_root, "/mean_matrices")

if (dir.exists(out_dir)) {
  print(paste(out_dir, "already exists"))
} else {
  dir.create(out_dir, recursive = TRUE)
  print(paste(out_dir, "created"))
}

################### 
# Define functions 
###################
# function to load RData connectivity matrices
read_rdata_matrix <- function(path, atlas) {
  m <- readRDS(path)
  matrix_subset_name <- grep(atlas, names(m), value = TRUE)
  m_atlas <- m[[matrix_subset_name]]
  if (is.data.frame(m_atlas)) m_atlas <- as.matrix(m_atlas)
  if (!is.matrix(m_atlas)) stop(path, " is not a matrix")
  storage.mode(m_atlas) <- "double"
  m_atlas
}

# Function to compute group mean
compute_mean_matrix <- function(subset_df, out_path, atlas) {
  if (nrow(subset_df) == 0) { message("No subjects for subset; skipping ", out_path); return(invisible(NULL)) }
  sum_mat <- NULL
  cnt_mat <- NULL
  
  for (p in subset_df$file) {
    m <- read_rdata_matrix(p, atlas)
    # initialize on first matrix; ensure all subsequent matrices have same dimensions
    if (is.null(sum_mat)) {
      sum_mat <- matrix(0, nrow = nrow(m), ncol = ncol(m))
      cnt_mat <- matrix(0L, nrow = nrow(m), ncol = ncol(m))
    } else if (!all(dim(m) == dim(sum_mat))) {
      stop("Matrix dimension mismatch at: ", p)
    }
    # create a logical mask not_na marking non-missing cells. set missing cells to 0 in m so they don’t affect the sum
    not_na <- !is.na(m)
    m[!not_na] <- 0
    
    # add the current matrix into the running sum and increase the count only at positions that were not NA
    sum_mat <- sum_mat + m
    cnt_mat[not_na] <- cnt_mat[not_na] + 1L 
  }
  
  mean_mat <- sum_mat # start with sum_mat as a template.
  zero <- cnt_mat == 0L # logical mask: TRUE wherever no subjects had a value 
  mean_mat[!zero] <- sum_mat[!zero] / cnt_mat[!zero] # for cells with at least one subject, compute means (sum / count)
  mean_mat[zero]  <- NA_real_ # wherever the count is zero (i.e., no subject had a value), set the mean to NA 
  
  saveRDS(mean_mat, out_path)
  message("Wrote ", out_path, " (", nrow(subset_df), " subjects)")
}

################### 
# Load files  
###################
# load demographics
demographics <- read.csv(config_data$demographics, stringsAsFactors = TRUE)
demographics <- demographics %>% select(sub, sex)

# get list of connectivity matrices
files <- list.files(mat_dir, pattern = "_ConnMatrices\\.RData$", full.names = TRUE)
if (length(files) == 0) stop("No *_ConnMatrices.RData found in: ", mat_dir)

# get subject ID = the filename prefix before "_ConnMatrices.RData"
get_id <- function(p) sub("_ConnMatrices\\.RData$", "", basename(p))
ids <- vapply(files, get_id, character(1), USE.NAMES = FALSE)

files_dt <- data.table(sub = ids, file = files)

df <- files_dt %>% inner_join(demographics, by = "sub")

if (nrow(df) == 0) stop("Join produced 0 rows—ID mismatch between filenames and demographics?")

######################### 
# Compute mean matrices
######################### 
female_df <- df[sex == "Female"]
male_df <- df[sex == "Male"]

compute_mean_matrix(female_df, file.path(out_dir, sprintf("%s_mean_connmatrix_female.rds", dataset)), atlas = atlas)
compute_mean_matrix(male_df, file.path(out_dir, sprintf("%s_mean_connmatrix_male.rds", dataset)), atlas = atlas)

print("Script finished!")


 
