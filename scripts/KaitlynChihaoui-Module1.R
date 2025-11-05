## Question 4
#Use any method to find either:
  
  #a simpler panel that achieves comparable classification accuracy

  #an alternative panel that achieves improved classification accuracy

#Benchmark your results against the in-class analysis.

# Lasso Method to find Simpler Panel w Comparable Accuracy

# Convert to matrices for glmnet
X_train <- training(biomarker_split) %>%
  select(-class) %>%
  as.matrix()

y_train <- training(biomarker_split) %>%
  pull(class) %>%
  as.numeric()

X_test <- testing(biomarker_split) %>%
  select(-class) %>%
  as.matrix()

y_test <- testing(biomarker_split) %>%
  pull(class)

# Fit LASSO with cross-validation
set.seed(101422)
cv_lasso <- cv.glmnet(X_train, y_train,
                      family = "binomial",
                      alpha = 1,
                      nfolds = 10)

# Extract coefficients at lambda.min
lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_proteins <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_proteins <- lasso_proteins[lasso_proteins != "(Intercept)"]


lasso_selected_coef <- lasso_coef[lasso_proteins, 1]
print(sort(abs(lasso_selected_coef), decreasing = TRUE))

# Make predictions on test set
lasso_pred_prob <- predict(cv_lasso,
                           newx = X_test,
                           s = "lambda.min",
                           type = "response")[,1]

# Evaluate using same style as in-class
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

lasso_test_results <- tibble(
  class = factor(y_test, labels = c("TD", "ASD")),
  pred = lasso_pred_prob,
  pred.group = factor(lasso_pred_prob > 0.5, labels = c("TD", "ASD"))
) %>%
  class_metrics(
    truth = class,
    estimate = pred.group,
    pred,
    event_level = "second"
  )

lasso_test_results


# Alternative Panel using Ridge Regression for Higher Accuracy
library(glmnet)
library(tidyverse)
library(tidymodels)

# Convert to matrices for glmnet
X_train <- training(biomarker_split) %>%
  select(-class) %>%
  as.matrix()

y_train <- training(biomarker_split) %>%
  pull(class) %>%
  as.numeric()

X_test <- testing(biomarker_split) %>%
  select(-class) %>%
  as.matrix()

y_test <- testing(biomarker_split) %>%
  pull(class)

# Fit Ridge with cross-validation
set.seed(101422)
cv_ridge <- cv.glmnet(X_train, y_train,
                      family = "binomial",
                      alpha = 0,
                      nfolds = 10)

# Extract coefficients at lambda.min
ridge_coef <- coef(cv_ridge, s = "lambda.min")
ridge_proteins <- rownames(ridge_coef)[which(ridge_coef != 0)]
ridge_proteins <- ridge_proteins[ridge_proteins != "(Intercept)"]


ridge_selected_coef <- ridge_coef[ridge_proteins, 1]
print(sort(abs(ridge_selected_coef), decreasing = TRUE))

# Make predictions on test set
ridge_pred_prob <- predict(cv_ridge,
                           newx = X_test,
                           s = "lambda.min",
                           type = "response")[,1]

# Evaluate using same style as in-class
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

ridge_test_results <- tibble(
  class = factor(y_test, labels = c("TD", "ASD")),
  pred = ridge_pred_prob,
  pred.group = factor(ridge_pred_prob > 0.5, labels = c("TD", "ASD"))
) %>%
  class_metrics(
    truth = class,
    estimate = pred.group,
    pred,
    event_level = "second"
  )

ridge_test_results

test_results
lasso_test_results
ridge_test_results

proteins_sstar
lasso_selected_coef
ridge_selected_coef

# Benchmark against results of in-class analysis

## Our in-class analysis yielded a classification accuracy of 0.839 using logistic regression, while our lasso regularization yielded a simpler panel of 6 
## proteins and a classification accuracy of 0.871, and our alternative method of ridge regression yielded a panel of 11 proteins at a classification 
## accuracy of 0.903. Lasso regularization gave us a comparable (within 3.2%) classification accuracy, while the ridge regression gave us a much higher (by 6.4%)
## classification accuracy.

## The selected proteins from the in-class analysis were: "DERM", "RELT", "Calcineurin", "IgD", "PTN", "FSTL1", "MAPK2", "TGF-b R III", "Notch 1", "ALCAM", 
## and "MATN2".  
## The selected proteins from the lasso panel were: DERM, Calcineurin, IgD, PTN, FSTL1, and MAPK2.
## The selected proteins from the ridge panel were: DERM, RELT, Calcineurin, IgD, PTN, FSTL1, MAPK2, TGF-b R III, Notch 1, ALCAM, and MATN2.

## Thus, the most selected proteins were the same between our in-class analysis and our ridge regression, but lasso only retained DERM, Calcineurin, IgD, 
## PTN, FSTL1, and MAPK2 from that subset.