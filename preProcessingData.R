# title: "R  Preprocessing and Download of Count Data - Breast Cancer TCGA"

# Installing and Loading Libraries            

if(!require("tidyverse")){install.packages("tidyverse")}

# Downloading TCGA-BRCA clinical data from Xenabrowser

# Survival data
# https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-BRCA.survival.tsv.gz"
destfile <- "Data/brca_survival.tsv.gz"
download.file(url, destfile)
brca.survi <- read_tsv(gzfile("Data/brca_survival.tsv.gz"))
brca.survi <- brca.survi %>% 
  mutate(sample = str_replace_all(sample, "-", ".")) %>% 
  column_to_rownames("sample") %>% 
  rename(status = OS, obs.time = OS.time, patient_id = '_PATIENT')

brca.survi <- as.data.frame(brca.survi)

# Clinical data
# https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-BRCA.GDC_phenotype.tsv.gz"
destfile <- "Data/brca_clinical.tsv.gz"
download.file(url, destfile)
brca.clini <- read_tsv(gzfile("Data/brca_clinical.tsv.gz"))
brca.clini <- brca.clini %>%
  dplyr::select(c("submitter_id.samples","prior_malignancy.diagnoses","age_at_initial_pathologic_diagnosis", "gender.demographic",
                  "sample_type_id.samples", "pathologic_M", "pathologic_N", "pathologic_T", "ethnicity.demographic", "race.demographic")) %>% 
  rename(sample = 'submitter_id.samples', 
         prior.dx = 'prior_malignancy.diagnoses', 
         age = 'age_at_initial_pathologic_diagnosis', 
         gender = 'gender.demographic',
         sample.type = 'sample_type_id.samples',
         metastasis = 'pathologic_M',
         neoplasm = 'pathologic_N',
         ajcc.stage = 'pathologic_T',
         ethnicity = 'ethnicity.demographic',
         race = 'race.demographic') %>% 
  mutate(sample.type = str_replace_all(sample.type, "01", "TP") ) %>% 
  mutate(sample.type = str_replace_all(sample.type, "11", "NT") ) %>% 
  filter(sample.type %in% c("TP", "NT")) %>%  
  mutate(sample = str_replace_all(sample, "-", ".")) %>% 
  filter(sample %in% row.names(brca.survi)) %>% 
  column_to_rownames(var = "sample") 

brca.clini <- cbind(brca.clini, brca.survi[rownames(brca.clini),])
brca.clini$codes <- rownames(brca.clini)

# Loading data form TGCA Ancestry - http://52.25.87.215/TCGAA/cancertype.php?cancertype=BRCA&pageSize=1000
brca.ances <- read_csv("Data/brca_tcga_ancestry.csv", na = "Not Available")
brca.ances <- brca.ances %>%
  rename(patient_id = "Patient ID", 
         cancer.type = "Cancer type", 
         race = 'Self-reported race', 
         ethnicity = 'Self-reported ethnicity',
         eigenstrat = 'EIGENSTRAT') %>% 
  column_to_rownames(var = "patient_id") 

# Combine the columns of brca.ances
brca.clini <- cbind(brca.clini, brca.ances[brca.clini$patient_id, ])

# Downloading TCGA-BRCA counts from Xenabrowser

# https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-BRCA.htseq_counts.tsv.gz"
destfile <- "brca_counts.tsv.gz"
download.file(url, destfile)
brca.count <- read_tsv(gzfile("C:/Users/raque/Downloads/Estagio/brca_counts.tsv.gz"))
brca.count <- as.data.frame(brca.count)
colnames(brca.count) <- gsub("-", "\\.", colnames(brca.count))
row.names(brca.count) <- sub("\\..*", "", brca.count$Ensembl_ID)
brca.count$Ensembl_ID <- NULL
brca.count <- 2^(brca.count)-1
brca.count <- round(brca.count, digits = 0)

# Select anotation dataset

# If you haven't already installed devtools...
# Use devtools to install the package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(annotables)

brca.annot <- grch38 %>%
  dplyr::filter(grch38$ensgene %in%  row.names(brca.count)) %>%
  dplyr::select(ensgene, symbol, description)

# Filtering Counts and Clinical data

brca.annot <- brca.annot[!duplicated(brca.annot$symbol), ]
brca.count <- brca.count[brca.annot$ensgene,]
rownames(brca.count) <- brca.annot$symbol
brca.clini <- brca.clini[rownames(brca.clini) %in% colnames(brca.count), ]
brca.count <- brca.count[,colnames(brca.count) %in% rownames(brca.clini)]

save(brca.count, brca.clini, brca.annot, file="C:/Users/raque/Downloads/Estagio/brca_count.RData", compress=T)
