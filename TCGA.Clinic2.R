######### TCGA Clinical Data Analysis #############
# checking variable relations to "over_surv_stt"
# using Tidyverse whenever possible


install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
library(tidyverse)
library(skimr)
library(finalfit)

setwd()


## 1. Importing data ---------------------

kirc_clin <- read_csv("kirc_clin.csv")


## 2. Taming data -------------------------
# use lubridate for dates
kirc_clin <- kirc_clin %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_diagnose = as.integer(year_diagnose))

glimpse(kirc_clin)
View(kirc_clin)


## 3. The dependent variable --------------------
# Check the number of levels. If greater than 2, run a ordinal logistic regression
# The independent variables are also called predicted or explanatory 
table(kirc_clin$over_surv_stt, exclude = NULL)


## 4. New predictors variables --------------------
# Making a categorical variable from a continuous one

hist(kirc_clin$age)
age <- kirc_clin[ , 'age']
age_group <- ifelse(age < 45, "under 45", 
                    ifelse(age >= 45 & age <= 65, "45 - 64",  
                           ifelse(age > 65 & age <= 75, "65 - 74",  
                                  ifelse(age > 75, "over 74", NA)))) 
# check if the new variable is correct
head(cbind(age_group, age)) 

kirc_clin$age_group <- as.factor(age_group) 
table(kirc_clin$age_group, exclude = NULL)


# Examine the candidate predictor variables.  

## 5. Categorical variables (chi-square test) --------------------
# Examine the relationship between chategorical variables and the dependent one, using cross tabulation.

# In the function prop.table(), using margin = 1 the table presents the row percentages, i.e. the frequencies of each column (over_surv_stt) in each row (explanatory variable). If we entered margin = 2, this would display the inverse, i.e. the frequencies of each row in each column.

# ToDo: table with all 
# t_metas <- table(kirc_clin$metastasis_stg, kirc_clin$over_surv_stt, exclude = NA)
# t_metas <- addmargins(round(100*prop.table(t_metas)))
# t_metas
# chisq.test(x = kirc_clin$metastasis_stg, y = kirc_clin$over_surv_stt) 


# summarise categorical variables by a categorical variable
explanatory_char <- kirc_clin %>% 
     select(-over_surv_stt) %>%
     select_if(is.factor) %>%
     names
dependent <- 'over_surv_stt'

table_char <- kirc_glm %>%
  summary_factorlist(dependent, explanatory_char, p=TRUE, 
                     na_include = FALSE, add_dependent_label=TRUE)
# To include NA into your statistic test: na_to_p = TRUE

table_char
# knitr::kable(table_char, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

warnings()
# Droping levels with narrow distributions -> check warnings ()
# Group some levels or drop one (NULL = 'level') when grouping is not possible  
# http://www.oncoguia.org.br/conteudo/estadiamento-do-cancer-de-rim/1811/240/

kirc_glm <- kirc_clin %>%
     mutate(neoplasm_stg = fct_collapse(neoplasm_stg, 'Stage I-II' = c('Stage I','Stage II'), 'Stage III-IV' = c('Stage III','Stage IV')),
            tumor_stg = fct_collapse(tumor_stg, 'T1-T2' = c('T1','T2'), 'T3-T4' = c('T3','T4')),
            histology_grd = fct_collapse(histology_grd, 'G1-G2' = c('G1','G2'), 'G3-G4' = c('G3','G4')),
            hemoglobin = fct_collapse(hemoglobin, 'Not.Low' = c('Normal', 'Elevated')),
            serum_ca = fct_collapse(serum_ca, 'Not.Low' = c('Normal', 'Elevated')),
            wbc = fct_collapse(wbc, 'Not.Elevated' = c('Low', 'Normal')))
                                                       
kirc_glm <- kirc_glm %>%
     mutate(histology_grd = fct_recode(histology_grd,  NULL = 'GX'),
            tumor_lateral = fct_recode(tumor_lateral, NULL = 'Bilateral'),
            race = fct_recode(race, NULL = 'Asian'))                                  
      

## 6. Saving dataset ------------------------------------
# Save dataset with different name if you processed any variable for glm analysis 
write_csv(kirc_glm, path = "your-path-here/kirc_glm.csv")

                            
## 7. Numeric variables (t.test) ---------------
# Examine the relationship between your candidate predictor variables
# For continuous variables, use pairwise correlations and scatterplot matrices

# kirc_clin_x <- kirc_clin %>%
#      group_by(over_surv_stt) %>%
#      select_if(is.numeric) %>%
#      summarise_at(vars(age:second_long_dim), mean, na.rm = TRUE) %>%
#      group_by(over_surv_stt)

      
# ggplot(kirc_clin, aes(age, fill= over_surv_stt)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(kirc_clin$age ~ kirc_clin$over_surv_stt) 
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=disease_free_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$disease_free_mth ~ kirc_clin$over_surv_stt) 
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=frac_genome_alter)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$frac_genome_alter ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(year_diagnose, fill= over_surv_stt)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(kirc_clin$year_diagnose ~ kirc_clin$over_surv_stt) 
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=long_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$long_dim ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=mutation_cnt)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$mutation_cnt ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=over_surv_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$over_surv_mth ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=short_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$short_dim ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=second_long_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$second_long_dim ~ kirc_clin$over_surv_stt)


# summarise categorical variables by a categorical variable
explanatory_num <- kirc_clin %>% 
     select_if(is.numeric) %>%
     names
dependent <- 'over_surv_stt'

table_fit_num <- kirc_clin %>%
     summary_factorlist(dependent, explanatory_num, p=TRUE, 
                        na_include = TRUE, add_dependent_label=TRUE)
# to include NA into your statistic test: na_to_p = TRUE

table_fit_num
# knitr::kable(table_fit, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

warnings()


# Correlation Matrix
# Use Pearson's (normal distribution) and Spearman (not-normal) correlations
# For continuous variables with no discontinuations and no obvious outliers
# Check collinearity (a strong linear relationship between predictors that causes problems in estimation of model parameters)
corr_num <- kirc_clin %>%
     select_if(is.numeric) %>%
     drop_na()

# Check the correlation between variables to exclude the higly correlated
cor_matrix <- cor(corr_num, method = "spearman")
cor_matrix <- round(cor_matrix, 2)
cor_matrix

# To visualize
pairs(~ age + disease_free_mth + frac_genome_alter + year_diagnose + long_dim + mutation_cnt + over_surv_mth + short_dim + second_long_dim, data=kirc_clin)

# Higlhy collinears: over_surv_mt & disease_free_mth 

