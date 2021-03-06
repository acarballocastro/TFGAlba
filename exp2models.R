library(tidyverse)
library(ggpubr)
library(shapr)
#library(ggforce)
library(caret)
library(ranger)
library(mgcv)
library(patchwork)

# Plots
source("extra/sina_plot.R")
source("extra/indiv_plot.R")

# Import and preprocessing

dataset <- read_csv("data/datasetADNI.csv")
head(dataset)

# Fixing covariables and response variable
x_var <- c("FDG","ABETA","PTAU","APOE4","PTGENDER","AGE","PTEDUCAT")
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

# Factor variables converted
# dataset$APOE4 = as.factor(dataset$APOE4)
# dataset$PTGENDER = as.factor(dataset$PTGENDER)
dataset$DXB = as.factor(dataset$DXB)
head(dataset)

# Train and test dataset
data_train <- dataset[train_index,]
data_test <- dataset[-train_index,]

# Model: Random forest ----

# NO da error pero tarda demasiado -> no tengo capacidad computacional

modelrf <- ranger(DXB ~ ., data = data_train, probability = TRUE)
print(modelrf)

# Model evaluation
# Prediction
pred = predict(modelrf, data_test)
pred_test = pred$predictions

# Converting to class
pred_test = colnames(pred_test)[max.col(pred_test, ties.method="first")]

# Creating confusion matrix
confusion = confusionMatrix(as.factor(c(y_test)), as.factor(pred_test))
print(confusion)

## Shapley values ----

### Symmetric ----

explainer_symmetric <- shapr(x_train, modelrf) 
p <- mean(y_train) # Expected prediction

#### Causal Shapley values ----
partial_order <- list(c(5,4,6,7), c(2), c(1,3))

explanation_causal <- explain(x_test, approach = "causal", 
                              explainer = explainer_symmetric, 
                              prediction_zero = p, ordering = partial_order,
                              confounding = c(TRUE, FALSE, TRUE), seed = 2022)

  


sina_causal <- sina_plot(explanation_causal)
# save limits of sina_causal plot for comparing against marginal and asymmetric
ylim_causal <- sina_causal$coordinates$limits$y

#### Marginal Shapley values ----
# Assumes one component with confounding
explanation_marginal <- explain(x_test, approach = "causal",
                                explainer = explainer_symmetric,
                                prediction_zero = p, ordering = list(c(1:7)),
                                confounding = TRUE, seed = 2020)

sina_marginal <- sina_plot(explanation_marginal) +
  coord_flip(ylim = ylim_causal) + ggtitle("Marginal Shapley values")

### Asymmetric ----

explainer_asymmetric <- shapr(x_train, modelrf, asymmetric = TRUE, ordering = partial_order)
p <- mean(y_train)

#### Asymmetric Shapley values ----

explanation_asymmetric <- explain(x_test, approach = "causal",
                                  explainer = explainer_asymmetric,
                                  prediction_zero = p, ordering = list(c(1:7)),
                                  confounding = TRUE, asymmetric = TRUE, seed = 2020)

sina_asymmetric <- sina_plot(explanation_asymmetric) +
  coord_flip(ylim = ylim_causal) + ggtitle("Asymmetric conditional Shapley values")

#### Asymmetric causal Shapley values ----

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

## Individual prediction explanations ----

source("extra/indiv_plot.R")

plot <- function(indiv, lim_inf = 0.25, lim_sup = 0.25){
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

plot(28)

# Model: Logistic regression for binary classification ----

# No puede tratar las variables como factores

dataset <- read_csv("data/datasetADNI.csv")
data_train <- dataset[train_index,]
data_test <- dataset[-train_index,]

# x and y with factors
x_train_df <- dataset[train_index, x_var]
y_train_df <- dataset[train_index, y_var]
x_test_df <- dataset[-train_index, x_var]
y_test_df <- dataset[-train_index, y_var]

modelglm <- glm(DXB ~ ., data = data_train, family = binomial)
print(modelglm)

# Model evaluation
# Prediction
pred_test = predict(modelglm, x_test_df, type = "response")

# Converting to class: c = 0.5
pred_test[(pred_test > 0.5)] = 1
pred_test[(pred_test < 0.5)] = 0

# Creating confusion matrix
confusion = confusionMatrix(as.factor(c(y_test)), as.factor(pred_test))
print(confusion)

## Shapley values ----

### Symmetric ----

explainer_symmetric <- shapr(x_train, modelglm) 
p <- mean(y_train) # Expected prediction

#### Causal Shapley values ----
partial_order <- list(c(5,4,6,7), c(2), c(1,3))

explanation_causal <- explain(x_test, approach = "causal", 
                              explainer = explainer_symmetric, 
                              prediction_zero = p, ordering = partial_order,
                              confounding = c(TRUE, FALSE, TRUE), seed = 2022)

sina_causal <- sina_plot(explanation_causal)
# save limits of sina_causal plot for comparing against marginal and asymmetric
ylim_causal <- sina_causal$coordinates$limits$y

#### Marginal Shapley values ----
# Assumes one component with confounding
explanation_marginal <- explain(x_test, approach = "causal",
                                explainer = explainer_symmetric,
                                prediction_zero = p, ordering = list(c(1:7)),
                                confounding = TRUE, seed = 2020)

sina_marginal <- sina_plot(explanation_marginal) +
  coord_flip(ylim = ylim_causal) + ggtitle("Marginal Shapley values")

### Asymmetric ----

explainer_asymmetric <- shapr(x_train, modelglm, asymmetric = TRUE, 
                              ordering = partial_order)
p <- mean(y_train)

#### Asymmetric Shapley values ----

explanation_asymmetric <- explain(x_test, approach = "causal",
                                  explainer = explainer_asymmetric,
                                  prediction_zero = p, ordering = list(c(1:7)),
                                  confounding = TRUE, asymmetric = TRUE, seed = 2020)

sina_asymmetric <- sina_plot(explanation_asymmetric) +
  coord_flip(ylim = ylim_causal) + ggtitle("Asymmetric conditional Shapley values")

#### Asymmetric causal Shapley values ----

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

## Individual prediction explanations ----

source("extra/indiv_plot.R")

plot <- function(indiv, lim_inf = 0.25, lim_sup = 0.25){
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

plot(28)
