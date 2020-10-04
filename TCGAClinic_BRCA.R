###### TCGA Clinical Data Analysis - BRCA #######
# cleaning and exploring dataset 

#### Loading Libraries ####
install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
library(tidyverse)
library(skimr)
library(finalfit)
library(dplyr)
library(ggplot2)

#### Download clinical data at combined_studies from cbioportal for Breast Cancer and TCGAA ####

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

# Visualizing
class(studies_brca_merge) 
dim(studies_brca_merge) 
names(studies_brca_merge) 
glimpse(studies_brca_merge)
skim(studies_brca_merge) 

#### Analyzing clinical data ####

# Frequencies of Each Status of Disease Status of Breast Cancer
studies_brca_merge %>% 
     count(DISEASE_FREE_STATUS) %>% 
     mutate(perc = n / nrow(studies_brca_merge)) -> subtype_brca

colnames(subtype_brca)  <- c("DISEASE_FREE_STATUS", "N", "Perc" )

suppressWarnings(
     ggplot(subtype_brca, aes(x = reorder(subtype_brca$DISEASE_FREE_STATUS, -subtype_brca$Perc) , y = subtype_brca$Perc )) + 
          scale_fill_brewer() + 
          geom_bar(aes(y = subtype_brca$Perc, fill = factor(subtype_brca$Perc)), colour = "black", stat = "identity") +
          geom_text(aes( label = scales::percent(subtype_brca$Perc),  fill = factor(..x..)), stat = "identity", vjust = -.4) +
          scale_y_continuous(labels = scales::percent) +
          geom_text(data = subtype_brca, aes(x = subtype_brca$DISEASE_FREE_STATUS, y = subtype_brca$Perc, label = subtype_brca$N), 
                    vjust = 1.8,  colour = "brown",  size = 4, show_guide = F) +
          guides(fill = FALSE) +
          ylab("Frequencies") +
          xlab("Disease Status of Breast Cancer"))

# Frequencies of Each Status of EIGENSTRAT of Breast Cancer

studies_brca_merge %>% 
     count(EIGENSTRAT) %>% 
     mutate(perc = n / nrow(studies_brca_merge)) -> subtype_brca

colnames(subtype_brca)  <- c("EIGENSTRAT", "N", "Perc" )

suppressWarnings(
     ggplot(subtype_brca, aes(x = reorder(subtype_brca$EIGENSTRAT, -subtype_brca$Perc) , y = subtype_brca$Perc )) + 
          scale_fill_brewer() + 
          geom_bar(aes(y = subtype_brca$Perc, fill = factor(subtype_brca$Perc)), colour = "black", stat = "identity") +
          geom_text(aes( label = scales::percent(subtype_brca$Perc),  fill = factor(..x..)), stat = "identity", vjust = -.4)+
          scale_y_continuous(labels = scales::percent) +
          geom_text(data = subtype_brca, aes(x = subtype_brca$EIGENSTRAT, y = subtype_brca$Perc, label = subtype_brca$N), 
                    vjust = 1.8,  colour = "brown",  size = 4, show_guide = F) +
          guides(fill = FALSE) +
          ylab("Frequencies") +
          xlab("Ancestry of Breast Cancer"))

# Frequencies of Each Status of SUBTYPE of Breast Cancer

studies_brca_merge %>% 
     count(ONCOTREE_CODE) %>% 
     mutate(perc = n / nrow(studies_brca_merge)) -> subtype_brca

colnames(subtype_brca)  <- c("ONCOTREE_CODE", "N", "Perc" )

suppressWarnings(
     ggplot(subtype_brca, aes(x = reorder(subtype_brca$ONCOTREE_CODE, -subtype_brca$Perc) , y = subtype_brca$Perc )) + 
          scale_fill_brewer() + 
          geom_bar(aes(y = subtype_brca$Perc, fill = factor(subtype_brca$Perc)), colour = "black", stat = "identity") +
          geom_text(aes( label = scales::percent(subtype_brca$Perc),  fill = factor(..x..)), stat = "identity", vjust = -.4) +
          scale_y_continuous(labels = scales::percent) +
          geom_text(data = subtype_brca, aes(x = subtype_brca$ONCOTREE_CODE, y = subtype_brca$Perc, label = subtype_brca$N), 
                    vjust = 1.8,  colour = "brown",  size = 4, show_guide = F) +
          guides(fill = FALSE) +
          ylab("Frequencies") +
          xlab("Oncotree of Breast Cancer"))

# Cumulative incidence of Disease Free/Recurred in BRCA Samples by Ancestry

library(survival)
library(survminer)

# Cumulative Incidence Function
require(cmprsk)

data <- data.frame(studies_brca_merge[,c("DISEASE_FREE_MONTHS_","DISEASE_FREE_STATUS","EIGENSTRAT")], check.names = TRUE, check.rows = TRUE)
fit_subtype <- cuminc(ftime = data$DISEASE_FREE_MONTHS, fstatus = data$DISEASE_FREE_STATUS, group = data$EIGENSTRAT)
fit_subtype

ggcompetingrisks(fit_subtype, palette = "Dark2",
                 legend = "top",
                 ggtheme = theme_bw())

# Downloading BRCA Samples by Normal Tissue and Tumour 

query.brca <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       platform = "Illumina HiSeq", 
                       file.type  = "normalized_results",
                       experimental.strategy = "RNA-Seq",
                       legacy = TRUE)

GDCdownload(query.brca, method = "api", files.per.chunk = 10)
data.brca <- GDCprepare(query.brca)

listPatient_EA <- as.vector(studies_brca_merge$PATIENT_ID[studies_brca_merge$EIGENSTRAT == "EA"])

listPatient_AA <- as.vector(studies_brca_merge$PATIENT_ID[studies_brca_merge$EIGENSTRAT == "AA"])

samples_EA <- query.brca$results[[1]]$cases[substr(query.brca$results[[1]]$cases,1,12) %in%  listPatient_EA]

samples_AA <- query.brca$results[[1]]$cases[substr(query.brca$results[[1]]$cases,1,12) %in%  listPatient_AA]

# Tumor 
samples_EA_TP <- TCGAquery_SampleTypes(barcode = samples_EA, typesample = "TP")
samples_AA_TP <- TCGAquery_SampleTypes(barcode = samples_AA, typesample = "TP")
# Normal Tissue
samples_EA_NT <- TCGAquery_SampleTypes(barcode = samples_EA, typesample = "NT")
samples_AA_NT <- TCGAquery_SampleTypes(barcode = samples_AA, typesample = "NT")