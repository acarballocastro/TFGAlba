library(tidyverse)
library(ggpubr)
library(shapr)
library(caret)
library(randomForest)

# Installing the causally enhanced shapr package from source
devtools::install_local("shapr-master", dependencies = TRUE)

dataset <- read_csv("data/datasetADNI.csv")
dataset <- dataset[,-1] # Removing ID
head(dataset)
dataset$APOE4 = as.factor(dataset$APOE4)
dataset$PTGENDER = as.factor(dataset$PTGENDER)
dataset$DX = as.factor(dataset$DX)
head(dataset)
dataset$PTGENDER = unclass(dataset$PTGENDER)-1
dataset$DX = unclass(dataset$DX)
head(dataset)

# Fixing covariables and response variable
x_var <- c("FDG","ABETA","PTAU","APOE4","PTGENDER","AGE","PTEDUCAT")
y_var <- "DX"

# Splitting in train-test (80%-20%) ----
set.seed(2022)
train_index <- caret::createDataPartition(dataset$DX, p = .8, list = FALSE, times = 1)

# Training data
x_train <- as.matrix(dataset[train_index, x_var])
y_train_nc <- as.matrix(dataset[train_index, y_var]) # not centered
y_train <- y_train_nc - mean(y_train_nc) # Tengo covariables categóricas

# Test data
x_test <- as.matrix(dataset[-train_index, x_var])
y_test_nc <- as.matrix(dataset[-train_index, y_var]) # not centered
y_test <- y_test_nc - mean(y_train_nc) # Tengo covariables categóricas

# Random forest ----
model <- randomForest(x = x_train, y = y_train_nc, ntree=500) 
print(model)

# Shapley values ----

explainer_symmetric <- shapr(x_train, model)                    
p <- mean(y_train)

# a. We compute the causal Shapley values on a given partial order (see paper)
partial_order <- list(1, c(2, 3), c(4:7)) # Poner el nuevo

explanation_causal <- explain(
  x_test,
  approach = "causal",
  explainer = explainer_symmetric,
  prediction_zero = p,
  ordering = partial_order,
  confounding = c(FALSE, TRUE, FALSE),
  seed = 2020
)

sina_causal <- sina_plot(explanation_causal)
# save limits of sina_causal plot for comparing against marginal and asymmetric
ylim_causal <- sina_causal$coordinates$limits$y

# b. For computing marginal Shapley values, we assume one component with confounding
explanation_marginal <- explain(
  x_test,
  approach = "causal",
  explainer = explainer_symmetric,
  prediction_zero = p,
  ordering = list(c(1:7)),
  confounding = TRUE,
  seed = 2020
)

sina_marginal <- sina_plot(explanation_marginal) +
  coord_flip(ylim = ylim_causal) + ylab("Marginal Shapley value (impact on model output)")

# c. Finally, we compute the asymmetric Shapley values for the same partial order
explainer_asymmetric <- shapr(x_train, model, asymmetric = TRUE, ordering = partial_order)
p <- mean(y_train)

explanation_asymmetric <- explain(
  x_test,
  approach = "gaussian",
  explainer = explainer_asymmetric,
  prediction_zero = p,
  ordering = partial_order,
  asymmetric = TRUE,
  seed = 2020
)

sina_asymmetric <- sina_plot(explanation_asymmetric) +
  coord_flip(ylim = ylim_causal) + ylab("Asymmetric conditional Shapley value (impact on model output)")

# d. Asymmetric causal Shapley values (very similar to the conditional ones)

explanation_asymmetric_causal <- explain(
  x_test,
  approach = "causal",
  explainer = explainer_asymmetric,
  prediction_zero = p,
  asymmetric = TRUE,
  ordering = partial_order,
  confounding = c(FALSE, TRUE, FALSE),
  seed = 2020
)

sina_asymmetric_causal <- sina_plot(explanation_asymmetric_causal) +
  coord_flip(ylim = ylim_causal) + ylab("Asymmetric causal Shapley value (impact on model output)")

# Faltan barplots y shapley value scatter plots
