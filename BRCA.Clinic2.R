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

brca_clin <- read_csv("~/Google Drive/BRCA-TCGAA/Data/brca_clin.csv")


## 2. Taming data -------------------------
# use lubridate for dates
brca_clin <- brca_clin %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_diagnose = as.integer(year_diagnose))

glimpse(brca_clin)
View(brca_clin)


## 3. The dependent variable --------------------
# Check the number of levels. If greater than 2, run a ordinal logistic regression
# The independent variables are also called predicted or explanatory 
table(brca_clin$over_surv_stt, exclude = NULL)
table(brca_clin$race, exclude = NULL)


## 4. Categorical variables (chi-square test) --------------------
# Examine the relationship between chategorical variables and the dependent one, using cross tabulation.

# In the function prop.table(), using margin = 1 the table presents the row percentages, i.e. the frequencies of each column (over_surv_stt) in each row (explanatory variable). If we entered margin = 2, this would display the inverse, i.e. the frequencies of each row in each column.

# ToDo: table with all 
t_race <- table(brca_clin$race, brca_clin$over_surv_stt, exclude = NA)
t_race <- addmargins(round(100*prop.table(t_race)))
t_race
chisq.test(x = brca_clin$race, y = brca_clin$over_surv_stt)


# summarise categorical variables by a categorical variable
explanatory_char <- brca_clin %>% 
     select_if(is.factor) %>%
     select(-race) %>%
     names
dependent <- 'race'

table_char <- brca_clin %>%
  summary_factorlist(dependent, explanatory_char, p=TRUE, 
                     na_include = FALSE, add_dependent_label=TRUE)
# To include NA into your statistic test: na_to_p = TRUE

table_char
# knitr::kable(table_char, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

warnings()
# Droping levels with narrow distributions -> check warnings ()
# Group some levels or drop one (NULL = 'level') when grouping is not possible  

                            
## 6. Numeric variables (t.test) ---------------
# Examine the relationship between your candidate predictor variables
# For continuous variables, use pairwise correlations and scatterplot matrices

# brca_clin_x <- brca_clin %>%
#      group_by(over_surv_stt) %>%
#      select_if(is.numeric) %>%
#      summarise_at(vars(age:over_surv_mth), mean, na.rm = TRUE) %>%
#      group_by(over_surv_stt)

# ggplot(brca_clin, aes(age, fill= over_surv_stt)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(brca_clin$age ~ brca_clin$over_surv_stt)

# ggplot(brca_clin, aes(x=over_surv_stt, y=disease_free_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(brca_clin$disease_free_mth ~ brca_clin$over_surv_stt)


# summarise numeric variables by a categorical variable
explanatory_num <- brca_clin %>% 
     select_if(is.numeric) %>%
     select(-year_diagnose) %>%
     names
dependent <- 'race'

table_fit_num <- brca_clin %>%
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
corr_num <- brca_clin %>%
     select_if(is.numeric) %>%
     drop_na()

# Check the correlation between variables to exclude the higly correlated
cor_matrix <- cor(corr_num, method = "spearman")
cor_matrix <- round(cor_matrix, 2)
cor_matrix

# To visualize
pairs(~ age + disease_free_mth + frac_genome_alter + year_diagnose  + mutation_cnt + over_surv_mth, data=brca_clin)
# Higlhy collinears: over_surv_mt & disease_free_mth 


## 7. Saving dataset ------------------------------------
# Save dataset with different name if you processed any variable for glm analysis 
write_csv(brca_glm, path = "your-path-here/brca_glm.csv")