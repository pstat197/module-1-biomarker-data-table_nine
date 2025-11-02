# repeat the analysis but carry out the entire selection procedure on a training partition -- 
# in other words, set aside some testing data at the very beginning and don't use it until you are 
# evaluating accuracy at the very end

# Necessary packages and seed and loading dataset
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
setwd("/Users/aarti/Downloads/Capstone/197a/module1/module-1-biomarker-data-table_nine")
load("data/biomarker-clean.RData")
set.seed(101422)

# Splitting the dataset
data_split <- initial_split(biomarker_clean, prop = 0.8, strata = group)
bio_train <- training(data_split)
bio_test <- testing(data_split)

# Selection procedure on the training set
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

ttests_out <- bio_train %>%
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
predictors <- bio_train %>%
  select(-c(group, ados))

response <- bio_train %>% pull(group) %>% factor()

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

biomarker_sstar <- bio_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar, 
           family = 'binomial')

# evaluate errors on test set
bio_test_sstar <- bio_test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = if_else(group == "ASD", "ASD", "TD")) %>%  # use character labels
  select(-group)

class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

# Get predictions
bio_test_pred <- bio_test_sstar %>%
  add_predictions(fit, type = "response") %>%
  rename(prob = pred) %>%   # rename the prediction column
  mutate(
    pred_class = if_else(prob > 0.5, "ASD", "TD"),
    class = factor(class, levels = c("TD", "ASD")),
    pred_class = factor(pred_class, levels = c("TD", "ASD"))
  )

# Evaluate performance
results <- bio_test_pred %>%
  class_metrics(
    truth = class,
    estimate = pred_class,
    prob,
    event_level = "second"
  )

print(results)

# The results are showing the results of specificity and accuracy and the area under the curve is 71.5% which 
# is pretty good. 
