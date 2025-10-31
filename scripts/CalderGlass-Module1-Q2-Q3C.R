library(tidyverse)
# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_dirty <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())

# export as r binary
save(list = 'biomarker_dirty', 
     file = 'data/biomarker-notrim.RData')

notrim = get(load('data/biomarker-notrim.RData'))

# 2.

# Constructed a second version of the dataset so that only the numeric columns
# would be present

notrim2 = subset(notrim, select = -c(group, ados))

# The scaling that was present in the original preprocessing has been done
# earlier, so this function just checks if the absolute value of each element
# in a given protein column is greater than 3 and returns the indices where
# it was true. 

outlier_detection = function(df){
  result = which(abs(df) > 3)
  return(result)
}

# Lapply allows for the function to be applied to every protein column in the
# dataset

indices = lapply(notrim2, outlier_detection)

# Get every index/subject where at least one of their protein levels was an
# outlier

unique_indices = unique(unlist(indices))
print(unique_indices)

# These are nearly all the indices with the exception of indices 45 and 110
# that had at least one outlying value for a given protein. 
# So the proportions of ASD vs TD for the outlying values is 75/152 ASD vs 
# 77/152 TD (so nearly 1:1)

# Setting up the indices so the counts for each index can be taken

not_unique_indices = unlist(indices)

# Indices that appear multiple times in this list would imply that they had
# outliers for more than one protein.

# Builds the table which shows how frequently each index had a protein level 
# that was an outlier.

counts_outliers = table(not_unique_indices)

View(counts_outliers)

which(counts_outliers >= 100)

# Each element in the bottom row of numbers is the index in the table
# While each element in the top row of numbers is the actual index of the subject
# These subjects have outlier protein levels in at least 100 of the 1317 proteins
# Of these 6 subjects with the most outliers in terms of protein levels,
# subjects 9 and 52 are of the Autism Spectrum Disorder (ASD) group while
# the remaining 4 subjects are of the TD group.

# 3 (b).

library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

proteins_s1

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

proteins_s2

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# Created test set
biomarker_test = testing(biomarker_split)


# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

test_results = biomarker_test %>%
  add_predictions(fit, type = 'response') %>% 
  mutate(class = factor(class, labels = c("TD", "ASD")),
         pred.group = factor(pred > 0.5, labels = c("TD", "ASD")))

test_results %>% 
  class_metrics(
    truth = class,
    estimate = pred.group,
    pred,
    event_level = "second"
  )

# When the number of significant proteins to be selected as a result of the 
# multiple t-testing increased from 10 to 20 proteins, the additional proteins
# that were selected were MAPK2, TGF-b R III, DAF, MMP-2, gp130 soluble, Notch
# 1, MIA, ALCAM, MATN2, and ROR1.

# The confusion matrix for the classification errors for ASD and TD remained
# the same for multiple t-testing after including 10 additional proteins.

# When the number of significant proteins to be selected as a result of the 
# random forest method increased from 10 to 20 proteins, the additional proteins
# that were selected were MAPK2, CK-MB, RET, Calcineurin, TSP4, Notch 1, PTN, 
# ERBB1, MATN2, and CSK.

# The sensitivity estimate decreased from 0.875 to 0.812. Accuracy stayed the
# same while specificity and roc_auc increased from 0.8 to 0.867 and 0.908 to 
# 0.946 respectively.


