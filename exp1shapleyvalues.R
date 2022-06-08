library(tidyverse)
library(ggpubr)
library(shapr)
#library(ggforce)
library(caret)
library(xgboost)
library(patchwork)

# Installing the causally enhanced shapr package from source
# devtools::install_local("shapr-master", dependencies = TRUE)

# Plots
source("extra/sina_plot.R")
source("extra/indiv_plot.R")

# Import and preprocessing ----

dataset <- read_csv("data/datasetADNI.csv")
head(dataset)

# Fixing covariables and response variable
x_var <- c("FDG","ABETA","PTAU","APOE4", "PTGENDER","AGE","PTEDUCAT")
y_var <- "DXB" # Binary classification

# Splitting in train-test (80%-20%)
set.seed(2022)
train_index <- caret::createDataPartition(dataset$DXB, p = .8, list = FALSE, times = 1)

# Training data
x_train <- as.matrix(dataset[train_index, x_var])
y_train <- as.matrix(dataset[train_index, y_var])

# Test data
x_test <- as.matrix(dataset[-train_index, x_var])
y_test <- as.matrix(dataset[-train_index, y_var]) 

# Model: XGBoost
# XGBoost requires classes to be in integer format

# Hyperparameter tuning
# Parameters: binary/logistic classification (supported by shapr)
# https://cran.r-project.org/web/packages/ParBayesianOptimization/vignettes/tuningHyperparameters.html

# hyper_tuning <- function(max_depth, min_child_weight, subsample) {
#   
# params = list(
#   objective="binary:logistic",
#   eval_metric="auc",
#   max_depth = max_depth, 
#   min_child_weight = min_child_weight, 
#   subsample = subsample
# )
# 
# # 10-fold cross validation
# xgbcv <- xgb.cv(params = params, data = x_train, label = y_train, 
#                 nrounds = 100, nfold = 10, prediction = T, showsd = T, 
#                 early.stop.round = 20, maximize = F, verbose = 0)
# 
# return(
#   list( 
#     Score = max(xgbcv$evaluation_log$test_auc_mean)
#     , nrounds = xgbcv$best_iteration
#   )
# )
# }
# 
# library("ParBayesianOptimization")
# bounds <- list( 
#   max_depth = c(2L, 10L)
#   , min_child_weight = c(1, 25)
#   , subsample = c(0.25, 1)
# )
# 
# optObj <- bayesOpt(
#   FUN = hyper_tuning
#   , bounds = bounds
#   , initPoints = 4
#   , iters.n = 3
# )
# 
# optObj$scoreSummary
# getBestPars(optObj)

# Model training

params = list(
  objective="binary:logistic",
  eval_metric="error"
)

modelxgb <- xgboost(data = x_train, label = y_train, nround = 100, 
                    verbose = FALSE, params = params)
print(modelxgb)

# Model evaluation
# Prediction
pred_test = predict(modelxgb, x_test)

# Converting to class: c = 0.5
pred_test[(pred_test >= 0.5)] = 1
pred_test[(pred_test < 0.5)] = 0

# Creating confusion matrix
confusion = confusionMatrix(as.factor(c(y_test)), as.factor(pred_test))
print(confusion)

# Shapley values ----

## Dependent features ----
explainer <- shapr(x_train, modelxgb) 
p <- mean(y_train) # Expected prediction

### Gaussian ----
explanation_gaussian <- explain(x_test, approach = "gaussian", 
                              explainer = explainer, prediction_zero = p, 
                              seed = 2022)

# Removing categorical variables
gaussian_cuant = NULL
gaussian_cuant$x_test = explanation_gaussian$x_test[,-c(4,5)]
gaussian_cuant$dt = explanation_gaussian$dt[,-c(5,6)]
gaussian_cuant$p = explanation_gaussian$p

sina_gaussian <- sina_plot(gaussian_cuant) + 
  ggtitle("Shapley values\nGaussian approach")
# save limits of sina_gaussian plot for comparing against marginal and asymmetric
ylim_gaussian <- sina_gaussian$coordinates$limits$y

### Copula ----
explanation_copula <- explain(x_test, approach = "copula", 
                                explainer = explainer, prediction_zero = p, 
                                seed = 2022)

# Removing categorical variables
copula_cuant = NULL
copula_cuant$x_test = explanation_copula$x_test[,-c(4,5)]
copula_cuant$dt = explanation_copula$dt[,-c(5,6)]
copula_cuant$p = explanation_copula$p

sina_copula <- sina_plot(copula_cuant) +
  coord_flip(ylim = ylim_gaussian) + ggtitle("Shapley values\nCopula approach")

### Empirical ----
explanation_empirical <- explain(x_test, approach = "empirical", 
                              explainer = explainer, prediction_zero = p, 
                              seed = 2022)

sina_empirical <- sina_plot(explanation_empirical) +
  coord_flip(ylim = ylim_gaussian) + ggtitle("Shapley values\nEmpirical approach")

## Symmetric ----

explainer_symmetric <- shapr(x_train, modelxgb) 
p <- mean(y_train) # Expected prediction

### Causal Shapley values ----
partial_order <- list(c(5,4,6,7), c(2), c(1,3))

explanation_causal <- explain(x_test, approach = "causal", 
                              explainer = explainer_symmetric, 
                              prediction_zero = p, ordering = partial_order,
                              confounding = c(TRUE, FALSE, TRUE), seed = 2022)

sina_causal <- sina_plot(explanation_causal)
# save limits of sina_causal plot for comparing against marginal and asymmetric
ylim_causal <- sina_causal$coordinates$limits$y

### Marginal Shapley values ----
# Assumes one component with confounding
explanation_marginal <- explain(x_test, approach = "causal",
                                explainer = explainer_symmetric,
                                prediction_zero = p, ordering = list(c(1:7)),
                                confounding = TRUE, seed = 2020)

sina_marginal <- sina_plot(explanation_marginal) +
  coord_flip(ylim = ylim_causal) + ggtitle("Marginal Shapley values")

## Asymmetric ----

explainer_asymmetric <- shapr(x_train, modelxgb, asymmetric = TRUE, ordering = partial_order)
p <- mean(y_train)

### Asymmetric Shapley values ----

explanation_asymmetric <- explain(x_test, approach = "causal",
                                  explainer = explainer_asymmetric,
                                  prediction_zero = p, ordering = list(c(1:7)),
                                  confounding = TRUE, asymmetric = TRUE, seed = 2020)

sina_asymmetric <- sina_plot(explanation_asymmetric) +
  coord_flip(ylim = ylim_causal) + ggtitle("Asymmetric conditional Shapley values")

### Asymmetric causal Shapley values ----

explanation_asymmetric_causal <- explain(x_test, approach = "causal",
                                         explainer = explainer_asymmetric,
                                         prediction_zero = p, asymmetric = TRUE,
                                         ordering = partial_order,
                                         confounding = c(TRUE, FALSE, TRUE),
                                         seed = 2020)

sina_asymmetric_causal <- sina_plot(explanation_asymmetric_causal) +
  coord_flip(ylim = ylim_causal) + ggtitle("Asymmetric causal Shapley values")

# Showing all plots together
(sina_causal + sina_marginal) / (sina_asymmetric + sina_asymmetric_causal)

# Individual prediction explanations ----

plot <- function(indiv, lim_inf = 0.2, lim_sup = 0.35){
indiv_1 <- indiv_plot(explanation_causal, digits = 3, plot_phi0 = FALSE, 
                      index_x_test = c(indiv)) +
    coord_flip(ylim = c(-lim_inf, lim_sup))
indiv_2 <- indiv_plot(explanation_marginal, digits = 3, plot_phi0 = FALSE, 
                      index_x_test = c(indiv)) +
  coord_flip(ylim = c(-lim_inf, lim_sup)) + labs(title = "Marginal Shapley values")
indiv_3 <- indiv_plot(explanation_asymmetric, digits = 3, plot_phi0 = FALSE, 
                      index_x_test = c(indiv)) +
    coord_flip(ylim = c(-lim_inf, lim_sup)) + labs(title = "Asymmetric Shapley values")
indiv_4 <- indiv_plot(explanation_asymmetric_causal, digits = 3, 
           plot_phi0 = FALSE, index_x_test = c(indiv)) + 
  coord_flip(ylim = c(-lim_inf, lim_sup)) + labs(title = "Asymmetric causal Shapley values")

# Showing all plots together
(indiv_1 + indiv_2) / (indiv_3 + indiv_4)}

# predict(modelxgb, x_test)
plot(54)
