###### TCGA Clinical Data Analysis #######
# cleaning and exploring dataset 
# using Tidyverse whenever possible


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
kirc_clin_raw <- read_delim("kirc_tcga_clinical_data.tsv", "\t", 
                            escape_double = FALSE, 
                            trim_ws = TRUE)

class(kirc_clin_raw) 
dim(kirc_clin_raw) 
names(kirc_clin_raw) 
glimpse(kirc_clin_raw)
skim(kirc_clin_raw) 
#View(kirc_clin_raw)


## 2. Cleaning data ---------------------------
# Select variables based on NA count (> 50% complete is a good choice!)
# ToDo: function to select variables with % completeness  

NA_fifty <- dim(kirc_clin_raw)[1]/2

NA_sum <- colSums(is.na(kirc_clin_raw))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum %>%
     filter(NA_sum < NA_fifty)

kirc_clean <- kirc_clin_raw %>%
     select(one_of(NA_sum$variables))


## Remove duplicate observations
kirc_clean0 <- kirc_clean %>%
     distinct_at('Patient ID', .keep_all = TRUE)

## Remove numeric variables with unique observations  
# ToDo: function to select variables with unique observations? 
kirc_clean0 %>%
     select_if(is.numeric) %>%
     skim()

kirc_clean1 <-  kirc_clean0  %>%
     select(!c('Last Alive Less Initial Pathologic Diagnosis Date Calculated Day Value', 
               'Number of Samples Per Patient', 
               'Sample type id'))


# Remove character variables with unique observations 
kirc_clean1 %>%
     select_if(is.character) %>%
     skim()

kirc_clean2 <- kirc_clean1  %>%
     select(!c('Study ID', 'Cancer Type', 'Cancer Type Detailed', 
               'Neoplasm Histologic Type Name', 'ICD-10 Classification', 
               'International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code', 
               'Informed consent verified', 'Is FFPE', 'Oncotree Code', 'Sample Type', 'Tumor Tissue Site'))


# Remove character variables with similar information - check each one!
# kirc_clean2 %>%
#      select_if(is.character) %>%
#      skim()
table(kirc_clean2$`Overall Survival Status`, exclude = NULL)
table(kirc_clean2$`Patient's Vital Status`, exclude = NULL)

kirc_clean3 <- kirc_clean2  %>%
     select(!c('Sample ID', 'Other Patient ID', 'Other Sample ID', 'Pathology Report File Name', 'Pathology report uuid', "Patient's Vital Status"))


# Remove other variables not directly related to patient - check each one!
# kirc_clean3 %>%
#      select_if(is.character) %>%
#      skim()

kirc_clean4 <- kirc_clean3  %>%
     select(!c('Form completion date','International Classification of Diseases for Oncology, Third Edition ICD-O-3 Histology Code','Primary Lymph Node Presentation Assessment Ind-3', 'Tissue Prospective Collection Indicator', 'Tissue Retrospective Collection Indicator','Vial number'))


## 3. Changing variables names ---------------------------
# Using snake_style 

kirc_clean4 <- kirc_clean4 %>%
     rename(patient_id = 'Patient ID',
            age = 'Diagnosis Age',
            metastasis_stg = 'American Joint Committee on Cancer Metastasis Stage Code',
            lymph_stg = 'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code',
            neoplasm_stg = 'Neoplasm Disease Stage American Joint Committee on Cancer Code',
            tumor_stg = 'American Joint Committee on Cancer Tumor Stage Code',
            disease_free_mth = 'Disease Free (Months)',
            disease_free_stt = 'Disease Free Status',
            ethnicity = 'Ethnicity Category', 
            frac_genome_alter = 'Fraction Genome Altered',
            histology_grd = 'Neoplasm Histologic Grade',
            hemoglobin = 'Hemoglobin level',
            neoadj_therapy = 'Neoadjuvant Therapy Type Administered Prior To Resection Text',
            prior_cancer = 'Prior Cancer Diagnosis Occurence',
            year_diagnose = 'Year Cancer Initial Diagnosis',
            tumor_lateral = 'Primary Tumor Laterality',
            long_dim = 'Longest Dimension',
            mutation_cnt = 'Mutation Count',
            over_surv_mth = 'Overall Survival (Months)',
            over_surv_stt = 'Overall Survival Status',
            platelet = 'Platelet count',
            race = 'Race Category',
            serum_ca = 'Serum calcium level',
            gender = 'Sex',
            short_dim = 'Shortest Dimension',
            second_long_dim = 'Specimen Second Longest Dimension',
            tissue_site = 'Tissue Source Site',
            person_neoplasm_stt = 'Person Neoplasm Status',
            wbc = 'WBC')


## 4. Taming data ------------------------------------------
# Use lubridate for dates
kirc_clean4 <- kirc_clean4 %>%
     mutate_if(is.character, as.factor) %>%
     mutate(patient_id = as.character(patient_id),
            age = as.integer(age),
            year_diagnose = as.integer(year_diagnose))


## 5. Checking NA patterns -----------------------------
# Check distincts types of NAs: MCAR, MAR, MNAR
kirc_clean4  %>%
     missing_plot()

missing_glimpse(kirc_clean4)


## 6. Checking numeric variables -----------------------------
# Check data distribution, unplausible values, outliers.
# Never delete an unusual value if the value is a possible one. 
# Deleting unusual values will bias your results and cause you to underestimate the variability in the observations.

# ToDo: codigo para fazer plot de todas as variaveis numericas de uma vez. 

# Describing numeric variables with summary().
# If the median and mean are similar, the distribution is likely roughly symmetrical. 
# Otherwise, it will be skewed to the right or to the left.
kirc_clean4 %>%
     select_if(is.numeric) %>%
     summary()

# Histograms or density plots
ggplot(kirc_clean4, aes(age)) +
     geom_histogram(bins = 20, alpha = 0.8, color = "red")

ggplot(kirc_clean4, aes(year_diagnose)) +
     geom_density(color = "red")


# Boxplots 
# Inter quartil range (IQR) = Q3 — Q1
# whiskers = ±1.58 IQR / √n ∗ IQR, where ‘n’ = samples
# Outliers = values below or above min and max whiskers values, respectively
     
ggplot(kirc_clean4, aes(x ='', y=disease_free_mth)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$disease_free_mth)

ggplot(kirc_clean4, aes(x ='', y=frac_genome_alter)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$frac_genome_alter)

ggplot(kirc_clean4, aes(x ='', y=long_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$long_dim)

ggplot(kirc_clean4, aes(x ='', y=mutation_cnt)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$mutation_cnt)

ggplot(kirc_clean4, aes(x ='', y=over_surv_mth)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$over_surv_mth)

ggplot(kirc_clean4, aes(x ='', y=short_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$short_dim)

ggplot(kirc_clean4, aes(x ='', y=second_long_dim)) +
     geom_boxplot(width = .5) +
     geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc_clean4$second_long_dim)


## 7. Checking categorical variables --------------------------
# Check frequency, lables and levels 
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging

kirc_clean4 %>%
     select_if(is.factor) %>%
     summary() 

# To do: fct_drop to drop unused levels

# agregating levels
kirc_clin <- kirc_clean4 %>%
     mutate(tumor_stg = fct_collapse(tumor_stg,
                                     T1 = c('T1', 'T1a', 'T1b'),
                                     T2 = c('T2', 'T2a', 'T2b'),
                                     T3 = c('T3', 'T3a', 'T3b', 'T3c')))

kirc_clin <- kirc_clin %>%
     mutate(prior_cancer = fct_collapse(prior_cancer, 
                                        Yes = c('Yes', 'Yes, History of Prior Malignancy', 'Yes, History of Synchronous/Bilateral Malignancy')))

kirc_clin <- kirc_clin %>%
     mutate(gender = fct_collapse(gender, Male = c('MALE', 'Male')))

kirc_clin <- kirc_clin %>%
     mutate(tissue_site = fct_collapse(tissue_site,
                                       A = c('A3', 'AK', 'AS'),
                                       B = c('B0', 'B2', 'B4', 'B8', 'BP'),
                                       C = c('CJ', 'CW', 'CZ')),
            tissue_site = fct_lump(tissue_site, 3, other_level = "OTHERS"))

# changing level names
kirc_clin <- kirc_clin %>%
     mutate(ethnicity = fct_recode(ethnicity, 'hispanic/latino'='HISPANIC OR LATINO', 'not hispanic/latino'='NOT HISPANIC OR LATINO'),
            race = fct_recode(race, Asian='ASIAN', 'Black/African.american'='BLACK OR AFRICAN AMERICAN', White='WHITE'),
            person_neoplasm_stt = fct_recode(person_neoplasm_stt, Tumor.Free='TUMOR FREE', With.Tumor='WITH TUMOR'))

# kirc_clin <- kirc_clin %>%
#      mutate(gender = if_else(gender %in% c('Male', 'Female'), 1, 0))

kirc_clin %>%
     select_if(is.factor) %>%
     summary()

# to visualize
kirc_clin %>% 
count(tissue_site) %>% 
     knitr::kable()

kirc_clin %>% 
     ggplot(aes(x = tissue_site, fill = tissue_site)) + 
     geom_bar() +
     theme_bw(15) +
     xlab("tissue site") + 
     ylab("frequency") + 
     theme(legend.position = 'none')

kirc_clin %>% 
ggplot(aes(x = tissue_site, y = age, fill = tissue_site)) + 
     geom_boxplot() +
     theme_bw(15) +
     xlab("tissue site") + 
     ylab("age") + 
     facet_wrap(~ gender) +
     theme(legend.position = 'none')


## 8. Correcting and checking again -----------------------
kirc_clin$disease_free_mth[kirc_clin$disease_free_mth == -11.79] <- 11.79
kirc_clin$disease_free_mth[kirc_clin$disease_free_mth == -0.62] <- 0.62

skim(kirc_clin)


## 9. Saving dataset ------------------------------------
write_csv(kirc_clin, path = "your_path_here/kirc_clin.csv")

rm(kirc_clean3, kirc_clean2, kirc_clean1, kirc_clean0, kirc_clean, kirc_clin_raw, NA_sum, NA_fifty)
