

PNC <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample.csv")
PNC <- PNC %>%
  mutate(sex = case_when(
    sex == 1 ~ "Male",
    sex == 2 ~ "Female",
    TRUE ~ NA_character_  
  ))
PNC$sex <- factor(PNC$sex, levels = c("Female", "Male"), ordered = TRUE)
write.csv(PNC, "/cbica/projects/network_replication/covariate_analyses/sex_diff/input/PNC/sample_info/PNC_demographics_finalsample.csv", row.names=F)



HCPD <- read.csv("/cbica/projects/network_replication/input/HCPD/sample_selection/HCPD_demographics_finalsample.csv")
HCPD$age <- HCPD$interview_age/12
HCPD <- HCPD %>%
  mutate(sex = case_when(
    sex == "M" ~ "Male",
    sex == "F" ~ "Female",
    TRUE ~ NA_character_  
  ))
HCPD$sex <- factor(HCPD$sex, levels = c("Female", "Male"), ordered = TRUE)
write.csv(HCPD, "/cbica/projects/network_replication/covariate_analyses/sex_diff/input/HCPD/sample_info/HCPD_demographics_finalsample.csv", row.names=F)
 
 

NKI <- read.csv("/cbica/projects/network_replication/input/NKI/sample_selection/NKI_demographics_finalsample.csv")
NKI$sub <- gsub("sub-", "", NKI$sub)
NKI$sex <- factor(NKI$sex, levels = c("Female", "Male"), ordered = TRUE)
write.csv(NKI, "/cbica/projects/network_replication/covariate_analyses/sex_diff/input/NKI/sample_info/NKI_demographics_finalsample.csv", row.names=F)



HBN <- read.csv("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_demographics_finalsample.csv")
HBN$sex <- factor(HBN$sex, levels = c("Female", "Male"), ordered = TRUE)
write.csv(HBN, "/cbica/projects/network_replication/covariate_analyses/sex_diff/input/HBN/sample_info/HBN_demographics_finalsample.csv", row.names=F)

