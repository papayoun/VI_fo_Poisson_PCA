
rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(pacman) # For progress bar
library(torch)
source("utils_generating_data.R") # Creates experiment_params.rds
library(fastDummies)

# Data --------------------------------------------------------------------

Y_full <- read.table("data_sets/borneo/data_abundance.txt", 
                sep = ";", header = TRUE)
X_full <- read.table("data_sets/borneo/data_soil_characteristics_after_pca.txt", 
                     sep = ";", header = TRUE) %>% 
  arrange(Site) %>% 
  dummy_cols("Sol", remove_selected_columns = TRUE)
kept_sites <- pull(X_full, Site)
X <- select(X_full, - Site) %>% 
  as.matrix()
Y <- Y_full %>% 
  filter(Site %in% kept_sites) %>% 
  arrange(Site) %>% 
  select(-Site) %>% 
  as.matrix()

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_functions.R")

result_VI <- get_CAVI(Y = Y, 
                      X = X,
                      q = 7,
                      seed = 5,  
                      n_steps = 1000, 
                      batch_prop = .4,
                      debug = FALSE, 
                      amortize = TRUE,
                      amortize_in_Y = TRUE,
                      updates = c(Lambda = TRUE, Sigma = TRUE,
                                  Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                  Beta = TRUE, Z = TRUE),
                      get_ELBO_freq = 10)

save(result_VI, file = "result_VI_amortized_on_borneo.RData")

plot(result_VI$ELBOS[-c(1:800),])

apply(result_VI$params$Lambda$M, 1, var) %>% plot()
