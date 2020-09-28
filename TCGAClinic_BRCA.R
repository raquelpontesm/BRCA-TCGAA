###### TCGA Clinical Data Analysis - BRCA #######
# cleaning and exploring dataset 

#### Loading Libraries ####
install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
library(tidyverse)
library(skimr)
library(finalfit)

# Download clinical data at combined_studies from cbioportal for Breast Cancer and TCGAA

# setwd()
brca_clin_raw <- read.delim("combined_study_clinical_brca.tsv")
brca_ancestry <- read.delim("tcga_ancestry_brca.csv", sep = ",", na.strings = "Not Available")

# Changing variables names
colnames(brca_clin_raw) <- gsub("\\_\\_", "\\_", gsub("\\.", "\\_", toupper(colnames(brca_clin_raw))))
colnames(brca_ancestry) <- gsub("\\_\\_", "\\_", gsub("\\.", "\\_", toupper(colnames(brca_ancestry))))

# Selecting the important variables 
brca_clin_raw <- brca_clin_raw[, c("PATIENT_ID", "ONCOTREE_CODE", "CANCER_TYPE_DETAILED", "AGE_AT_DIAGNOSIS", "DISEASE_FREE_EVENT",  "DISEASE_FREE_MONTHS_", "OVERALL_SURVIVAL_MONTHS_", "OVERALL_SURVIVAL_STATUS", "DISEASE_FREE_STATUS", "TIME_TO_DEATH_MONTHS_", "TIME_TO_METASTASIS_MONTHS_", "METASTATIC_RECURRENCE_TIME", "PATIENT_S_VITAL_STATUS", "MEDR_TIME_TO_METASTATIC_DIAGNOSIS_CALCULATED_MONTHS_", "TYPE_OF_BREAST_SURGERY")]
brca_clin_raw <- brca_clin_raw[!is.na(brca_clin_raw$DISEASE_FREE_MONTHS_),]

# Combined studies
studies_brca_merge <- merge(brca_clin_raw, brca_ancestry, by = "PATIENT_ID")

# Remove duplicates
studies_brca_merge <- studies_brca_merge[!duplicated(studies_brca_merge$PATIENT_ID), ]

#Visualizing
class(studies_brca_merge) 
dim(studies_brca_merge) 
names(studies_brca_merge) 
glimpse(studies_brca_merge)
skim(studies_brca_merge) 

