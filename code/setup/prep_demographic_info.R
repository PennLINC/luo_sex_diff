
library(dplyr)
library(magrittr)

# Set project root depending on whether code is running on PARCC or locally
if (dir.exists("/ceph/projects/sattertt/pennlinc-parcc/network_replication")) {
  root <- "/ceph/projects/sattertt/pennlinc-parcc/network_replication"
} else if (dir.exists("~/parcc/network_replication")) {
  root <- path.expand("~/parcc/network_replication")
} else {
  stop("Could not find network_replication project.")
}

# Replication analysis directory
rep_root <- file.path(root, "covariate_analyses/sex_diff")


# PNC
PNC <- read.csv(file.path(root, "input/PNC/sample_selection/PNC_demographics_finalsample.csv"))

PNC <- PNC %>%
  mutate(sex = case_when(
    sex == 1 ~ "Male",
    sex == 2 ~ "Female",
    TRUE ~ NA_character_
  ))

PNC$sex <- factor(PNC$sex, levels = c("Male", "Female"))

PNC_out <- file.path(rep_root, "input/PNC/sample_info/PNC_demographics_finalsample.csv")
write.csv(PNC, PNC_out, row.names = FALSE)
cat("Wrote:", PNC_out, "\n")


# HCPD
HCPD <- read.csv(file.path(root, "input/HCPD/sample_selection/HCPD_demographics_finalsample.csv"))

HCPD$age <- HCPD$interview_age / 12

HCPD <- HCPD %>%
  mutate(sex = case_when(
    sex == "M" ~ "Male",
    sex == "F" ~ "Female",
    TRUE ~ NA_character_
  ))

HCPD$sex <- factor(HCPD$sex, levels = c("Male", "Female"))

HCPD_out <- file.path(rep_root, "input/HCPD/sample_info/HCPD_demographics_finalsample.csv")
write.csv(HCPD, HCPD_out, row.names = FALSE)
cat("Wrote:", HCPD_out, "\n")


# NKI
NKI <- read.csv(file.path(root, "input/NKI/sample_selection/NKI_demographics_finalsample.csv"))

NKI$sub <- gsub("sub-", "", NKI$sub)
NKI$sex <- factor(NKI$sex, levels = c("Male", "Female"))

NKI_out <- file.path(rep_root, "input/NKI/sample_info/NKI_demographics_finalsample.csv")
write.csv(NKI, NKI_out, row.names = FALSE)
cat("Wrote:", NKI_out, "\n")


# HBN
HBN <- read.csv(file.path(root, "input/HBN/sample_selection/HBN_demographics_finalsample.csv"))

HBN$sex <- factor(HBN$sex, levels = c("Male", "Female"))

HBN_out <- file.path(rep_root, "input/HBN/sample_info/HBN_demographics_finalsample.csv")
write.csv(HBN, HBN_out, row.names = FALSE)
cat("Wrote:", HBN_out, "\n")
 