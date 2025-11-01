# Preprocessing File:

library(tidyverse)
library(dplyr)


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
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
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
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

# export as r binary
save(list = 'biomarker_clean', 
     file = 'data/biomarker-clean.RData')






# Preliminary Analysis File:

library(tidyverse)
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
  slice_min(p.adj, n = 10) %>%
  pull(protein)

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
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

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

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  class_metrics(estimate = factor(pred > 0.5),
                truth = factor(class), pred,
                event_level = 'second')







# Question 1: 


library(tidyverse)

# Read in biomarker data
biomarker_raw <- read_csv("data/biomarker-raw.csv", skip = 1)

# Rename the first column for clarity
biomarker_raw <- biomarker_raw %>%
  rename(group = 1) %>%
  filter(!is.na(group))   # remove any rows without a group label

# Inspect a random sample of 5 protein columns
set.seed(102625)
sample_proteins <- sample(names(biomarker_raw)[-1], 5)

# Plot raw distributions for a sample of proteins
biomarker_raw %>%
  pivot_longer(cols = all_of(sample_proteins), 
               names_to = "protein", 
               values_to = "level") %>%
  ggplot(aes(x = level)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Distribution of Raw Protein Levels",
       x = "Raw Value",
       y = "Count")




# From the distribution of raw values for a sample of proteins, we can see that
# the data is heavily right-skewed. Therefore there are many small values of the
# proteins and very few larger values. Taking the log-transformation of the
# protein levels compresses the larger values and stretches the smaller values, 
# reducing right skew and making the data more symmetric and closer to a normal
# distribution. Additionally, using log-transformation can help to stabilize 
# the variance among protein levels. Thus, using the log-transformation on the 
# protein levels reduces skewness to make the distribution approximately normal
# and stabilizes variance, therefore helping the data to meet model assumptions 
# for tests including the t-test and logistic regression. 





# Question 3(c)

library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

# Load cleaned data
load('data/biomarker-clean.RData')

## MULTIPLE TESTING
####################
test_fn <- function(.df) {
  t_test(.df,
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = FALSE)
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = 'protein', values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m * hm * p_value / rank)

## RANDOM FOREST
##################
predictors <- biomarker_clean %>% select(-c(group, ados))
response <- biomarker_clean %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(x = predictors,
                       y = response,
                       ntree = 1000,
                       importance = TRUE)

## FUZZY INTERSECTION (replaces hard intersect)
###############################################

# Rank proteins by t-test (smallest p-value = best)
ranked_ttest <- ttests_out %>%
  select(protein, p_value) %>%
  mutate(rank_t = rank(p_value, ties.method = "average"))

# Rank proteins by Random Forest importance (largest Gini = best)
ranked_rf <- rf_out$importance %>%
  as_tibble(rownames = "protein") %>%
  select(protein, MeanDecreaseGini) %>%
  mutate(rank_rf = rank(-MeanDecreaseGini, ties.method = "average"))

# Combine rankings and compute average rank
combined_ranks <- full_join(ranked_ttest, ranked_rf, by = "protein") %>%
  mutate(avg_rank = (rank_t + rank_rf) / 2) %>%
  arrange(avg_rank)

# Select top proteins based on average rank (fuzzy intersection)
proteins_sstar <- combined_ranks %>%
  slice_min(avg_rank, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################
biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

fit <- glm(class ~ ., 
           data = training(biomarker_split),
           family = 'binomial')

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(
    pred_class = factor(pred > 0.5, levels = c(FALSE, TRUE)),
    truth = factor(class, levels = c(FALSE, TRUE))
  ) %>%
  class_metrics(
    estimate = pred_class,
    truth = truth,
    pred,
    event_level = "second"
  )



# Instead of using a hard intersection to combine sets of top predictive proteins
# across selection methods, we used a fuzzy intersection based on average rank. 
# We sorted the p-values of the t-tests from smallest to largest and the Gini 
# values of the random forest method from largest to smallest. Each protein
# received a rank according to these values. We then averaged these ranks and 
# selected the top 10 proteins with the smallest average ranks. This method
# captures proteins that are consistently strong across both selection methods,
# even if they don't appear in both top 10 sets.


# When using the fuzzy intersection method compared to a hard intersection, the
# sensitivity decreased from 0.812 to 0.75, the specificity increased from 0.733 
# to 0.867, the accuracy value increased from 0.774 to 0.806, and the roc_auc value 
# increased slightly from 0.883 to 0.892. Because the positive class corresponds 
# to ASD, this means the fuzzy intersection made the model slightly worse at 
# correctly identifying individuals with ASD (lower sensitivity), but better at 
# correctly identifying typically developing (TD) individuals (higher specificity). 
# Despite the drop in sensitivity, the overall accuracy and ROC AUC improved, 
# suggesting that including proteins that were highly ranked by at least one 
# selection method rather than both helped the model generalize slightly better 
# and improved its ability to discriminate between ASD and TD overall.