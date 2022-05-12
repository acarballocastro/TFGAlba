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
dataset <- dataset[,-1] # Removing ID
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

# Model: XGBoost
# XGBoost requires classes to be in integer format

# Parameters: binary/logistic classification (supported by shapr)
params = list(
  objective="binary:logistic",
  eval_metric="error"
)
modelxgb <- xgboost(data = x_train, label = y_train, nround = 100, 
                    verbose = FALSE, params = params)
print(modelxgb)
# predict(modelxgb, x_test)

# Shapley values ----

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

explanation_asymmetric <- explain(x_test, approach = "gaussian",
                                  explainer = explainer_asymmetric,
                                  prediction_zero = p, ordering = partial_order,
                                  asymmetric = TRUE, seed = 2020)

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

# Scatter plot ----

sv_correlation_df <- data.frame(
  valtemp = x_test[, "ABETA"],
  sv_marg_FDG = explanation_marginal$dt$FDG,
  sv_caus_FDG = explanation_causal$dt$FDG,
  sv_marg_ABETA = explanation_marginal$dt$ABETA,
  sv_caus_ABETA = explanation_causal$dt$ABETA
)

scatterplot_topleft <- 
  ggplot(sv_correlation_df, aes(x = sv_marg_ABETA, y = sv_marg_FDG, color = valtemp)) + 
  geom_point(size = 1)+xlab("MSV ABETA")+ylab("MSV FDG")+
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5))  + 
  scale_color_gradient(low="blue", high="red") +
  theme_minimal() + 
  theme(text = element_text(size = 12), 
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

scatterplot_topright <- 
  ggplot(sv_correlation_df, aes(x = sv_caus_FDG, y = sv_marg_ABETA, color = valtemp)) + 
  geom_point(size = 1) + scale_color_gradient(low="blue", high="red") +
  xlab("CSV ABETA") + ylab("MSV ABETA") + 
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) + 
  theme_minimal() +
  theme(text = element_text(size=12), axis.title.x = element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

scatterplot_bottomleft <- 
  ggplot(sv_correlation_df, aes(x = sv_marg_ABETA, y = sv_caus_ABETA, color = valtemp)) +
  geom_point(size = 1) + scale_color_gradient(low="blue", high="red") + 
  ylab( "CSV ABETA") + xlab("MSV ABETA") +  
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5))  + 
  theme_minimal() +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))

scatterplot_bottomright <- 
  ggplot(sv_correlation_df, aes(x = sv_caus_FDG, y = sv_caus_ABETA, color = valtemp)) +
  geom_point(size = 1) + ylab("CSV ABETA") + xlab( "CSV FDG") + 
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5))  + 
  scale_color_gradient(low="blue", high="red")+
  theme_minimal() +
  theme(text = element_text(size=12), axis.text.x=element_text(size=12),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

grid_top <- ggarrange(scatterplot_topleft, scatterplot_topright, legend = "none")
grid_bottom <- ggarrange(scatterplot_bottomleft, scatterplot_bottomright, legend = "none")

# Bar plot ----

# Get test set index for two data points with similar temperature
# 1. 2012-10-09 (October)
# 2. 2012-12-03 (December)

october <- which(as.integer(row.names(x_test)) == which(bike$dteday == "2012-10-09"))
december <- which(as.integer(row.names(x_test)) == which(bike$dteday == "2012-12-03"))

# predicted values for the two points
# predict(model, x_test)[c(october, december)] + mean(y_train_nc)

dt_marginal <- explanation_marginal$dt %>%
  dplyr::slice(c(october, december)) %>%
  select(cosyear, temp) %>%
  mutate(date = c("2012-10-09", "2012-12-03"), type = 'Marginal')

dt_causal <- explanation_causal$dt %>%
  dplyr::slice(c(october, december)) %>%
  select(cosyear, temp) %>%
  mutate(date = c("2012-10-09", "2012-12-03"), type = 'Causal')

dt_asymmetric <- explanation_asymmetric$dt %>%
  dplyr::slice(c(october, december)) %>%
  select(cosyear, temp) %>%
  mutate(date = c("2012-10-09", "2012-12-03"), type = 'Asymmetric')

dt_all <- dt_marginal %>% pivot_longer(c(cosyear, temp)) %>%
  rbind(dt_causal %>% pivot_longer(c(cosyear, temp))) %>%
  rbind(dt_asymmetric %>% pivot_longer(c(cosyear, temp)))

bar_plots <- ggplot(dt_all, aes(x = name, y = value, group = interaction(date, name), 
                                fill = date, label = round(value, 2))) +
  geom_col(position = "dodge") +
  theme_classic() + ylab("Shapley value") +
  facet_wrap(vars(type)) + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = c('indianred4', 'ivory4')) + 
  theme(legend.position = c(0.75, 0.25), axis.title = element_text(size = 20),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14))

