########## TCGA-BRCA Clinical Data ###########
# cleaning & exploring dataset with Tidyverse


install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
library(tidyverse)
library(skimr)
library(finalfit)

# setwd()


## 1. Data importing and visualizing ---------------------------
# Download clinical data at GDC (select primary site, project and choose 'clinical' at exploration page).
# Using read_csv or read_delim, na = c("", "NA")
brca_clin_raw <- read_delim("~/Google Drive/BRCA-TCGAA/Data/brca_tcga_clinical_data.tsv", "\t", 
                            escape_double = FALSE, 
                            trim_ws = TRUE)

class(brca_clin_raw) 
dim(brca_clin_raw) 
names(brca_clin_raw) 
glimpse(brca_clin_raw)
skim(brca_clin_raw) 
#View(brca_clin_raw)


## 2. Cleaning data ---------------------------
# Select variables based on NA count (> 50% complete is a good choice!)
# ToDo: function to select variables with % completeness  

NA_fifty <- dim(brca_clin_raw)[1]/2

NA_sum <- colSums(is.na(brca_clin_raw))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum %>%
     filter(NA_sum < NA_fifty)

brca_clean <- brca_clin_raw %>%
     select(one_of(NA_sum$variables))

## Remove duplicate observations
brca_clean1 <- brca_clean %>%
     distinct_at('Patient ID', .keep_all = TRUE)

## Numeric variables 
# Remove variables with unique observations and non-related ones
# ToDo: function to select variables with unique observations? 
brca_clean1 %>%
     select_if(is.numeric) %>%
     skim()

brca_clean2 <-  brca_clean1  %>%
     select(!c('Days to Sample Collection.',
               'Last Alive Less Initial Pathologic Diagnosis Date Calculated Day Value',
               'Positive Finding Lymph Node Hematoxylin and Eosin Staining Microscopy Count',
               'Lymph Node(s) Examined Number',
               'Number of Samples Per Patient',
               'Sample Initial Weight',
               'Sample type id'))

## Character variables 
# Remove character variables with unique observation, similar information, or not related - check each one! 
brca_clean2 %>%
     select_if(is.character) %>%
     skim()

brca_clean3 <- brca_clean2  %>%
     select(!c('Study ID','Sample ID','American Joint Committee on Cancer Publication Version Type',
               'Cancer Type', 'Cancer Type Detailed', 'Form completion date','HER2 ihc score',
               'Neoplasm Histologic Type Name','ICD-10 Classification', 
               'International Classification of Diseases for Oncology, Third Edition ICD-O-3 Histology Code',
               'International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code', 
               'Informed consent verified','Is FFPE',
               'Primary Lymph Node Presentation Assessment Ind-3',
               'First Pathologic Diagnosis Biospecimen Acquisition Method Type',
               'Other Patient ID', 'Other Sample ID', 'Pathology Report File Name',
               'Pathology report uuid','Disease Surgical Margin Status',
               'Patient Primary Tumor Site','Tissue Prospective Collection Indicator',
               'Tissue Retrospective Collection Indicator','Sample Type',
               'Sex','Tumor Tissue Site','Tissue Source Site','Vial number'))

skim(brca_clean3)

## Logical variables
brca_clean3 <- brca_clean3  %>%
     select(!'Oct embedded')


## 3. Changing variables names ---------------------------
# Using snake_style 

colnames(brca_clean3) <- tolower(colnames(brca_clean3))
colnames(brca_clean3) <- str_replace_all(colnames(brca_clean3), "[' ']", "_")

brca_clean3 <- brca_clean3 %>%
     rename(patient_id = 'patient_id',
            age = 'diagnosis_age',
            metastasis_stg = 'american_joint_committee_on_cancer_metastasis_stage_code',
            lymph_stg = 'neoplasm_disease_lymph_node_stage_american_joint_committee_on_cancer_code',
            neoplasm_stg = 'neoplasm_disease_stage_american_joint_committee_on_cancer_code',
            tumor_stg = 'american_joint_committee_on_cancer_tumor_stage_code',
            disease_free_mth = 'disease_free_(months)',
            disease_free_stt = 'disease_free_status',
            er_stt = 'er_status_by_ihc',
            ethnicity = 'ethnicity_category', 
            frac_genome_alter = 'fraction_genome_altered',
            neoadj_therapy = 'neoadjuvant_therapy_type_administered_prior_to_resection_text',
            prior_cancer = 'prior_cancer_diagnosis_occurence',
            her2_stt = 'ihc-her2',
            year_diagnose = 'year_cancer_initial_diagnosis',
            menopause = 'menopause_status',
            micromet = 'micromet_detection_by_ihc',
            mutation_cnt = 'mutation_count',
            oncotree = 'oncotree_code',
            over_surv_mth = 'overall_survival_(months)',
            over_surv_stt = 'overall_survival_status',
            pr_stt = 'pr_status_by_ihc',
            race = 'race_category',
            surgical = 'surgical_procedure_first',
            person_neo_stt = 'person_neoplasm_status')


## 4. Taming data ------------------------------------------
# Use lubridate for dates
brca_clean3 <- brca_clean3 %>%
     mutate_if(is.character, as.factor) %>%
     mutate(patient_id = as.character(patient_id))


## 5. Checking NA patterns -----------------------------
# Check distincts types of NAs: MCAR, MAR, MNAR
brca_clean3  %>%
     missing_plot()

missing_glimpse(brca_clean3)


## 6. Checking numeric variables -----------------------------
# Check data distribution, unplausible values, outliers.
# Never delete an unusual value if the value is a possible one. 
# Deleting unusual values will bias your results and cause you to underestimate the variability in the observations.
brca_clean3 %>%
     select_if(is.numeric) %>%
     summary()

brca_clean3$disease_free_mth[brca_clean3$disease_free_mth == -0.23] <- 0.23
brca_clean3$over_surv_mth[brca_clean3$over_surv_mth == -0.23] <- 0.23

# Graphics
ggplot(brca_clean3, aes(age)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

ggplot(brca_clean3, aes(disease_free_mth)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

ggplot(brca_clean3, aes(year_diagnose)) +
     geom_bar(stat="count")

ggplot(brca_clean3, aes(over_surv_mth)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

# Boxplots 
# Inter quartil range (IQR) = Q3 — Q1
# whiskers = ±1.58 IQR / √n ∗ IQR, where ‘n’ = samples
# Outliers = values below or above min and max whiskers values, respectively

ggplot(brca_clean3, aes(x ='', y=frac_genome_alter)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(brca_clean3$frac_genome_alter)

ggplot(brca_clean3, aes(x ='', y=mutation_cnt, na.rm = TRUE)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(brca_clean3$mutation_cnt)


## 7. Checking categorical variables --------------------------
# Check frequency, lables and levels 
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging
brca_clean3 %>%
select_if(is.factor) %>%
summary()

# agregating levels
brca_clean3 <- brca_clean3 %>%
     mutate(tumor_stg = fct_collapse(tumor_stg,
                                     T1 = c('T1','T1a','T1b','T1c'),
                                     T2 = c('T2','T2a','T2b'),
                                     T3 = c('T3','T3a'),
                                     T4 = c('T4','T4b','T4d')))


# changing level names
brca_clean3 <- brca_clean3 %>%
     mutate(disease_free_stt = fct_recode(disease_free_stt, DiseaseFree='0:DiseaseFree', RecurrProgress='1:Recurred/Progressed'),
          ethnicity = fct_recode(ethnicity, 'hispanic/latino'='HISPANIC OR LATINO', 'not hispanic/latino'='NOT HISPANIC OR LATINO'),
          menopause = fct_recode(menopause, Intermediate='Indeterminate (neither Pre or Postmenopausal)', Peri='Peri (6-12 months since last menstrual period)', Post='Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)', Pre='Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)'),
          over_surv_stt = fct_recode(over_surv_stt, living='0:LIVING', deceased='1:DECEASED'),
          race = fct_recode(race, Asian='ASIAN', 'Black/African.american'='BLACK OR AFRICAN AMERICAN', White='WHITE'),
          staging_system = fct_recode(staging_system, AxilLN='Axillary lymph node dissection alone', NoAxil='No axillary staging', Other='Other (specify)', SentinelAxil='Sentinel lymph node biopsy plus axillary dissection', SentinelLN='Sentinel node biopsy alone'),
          surgical = fct_recode(surgical, RadicalMastectomy='Modified Radical Mastectomy', SimpleMastectomy='Simple Mastectomy'),
          person_neo_stt = fct_recode(person_neo_stt, TumorFree='TUMOR FREE', WithTumor='WITH TUMOR'))


# to visualize
# brca_clean3 %>% 
# count(tissue_site) %>% 
#      knitr::kable()
# 
# brca_clean3 %>% 
#      ggplot(aes(x = tissue_site, fill = tissue_site)) + 
#      geom_bar() +
#      theme_bw(15) +
#      xlab("tissue site") + 
#      ylab("frequency") + 
#      theme(legend.position = 'none')
# 
# brca_clean3 %>% 
# ggplot(aes(x = tissue_site, y = age, fill = tissue_site)) + 
#      geom_boxplot() +
#      theme_bw(15) +
#      xlab("tissue site") + 
#      ylab("age") + 
#      facet_wrap(~ gender) +
#      theme(legend.position = 'none')

skim(brca_clean3)


## 8. Saving dataset ------------------------------------
write_csv(brca_clean3, path = "~/Google Drive/BRCA-TCGAA/Data/brca_clin.csv")

rm(brca_clean3, brca_clean2, brca_clean1, brca_clean, brca_clin_raw, NA_sum, NA_fifty)
