
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

result_VI <- NULL
if (file.exists("result_VI_full_CAVI_on_borneo.RData"))
  load("result_VI_full_CAVI_on_borneo.RData")
result_VI <- get_CAVI(Y = Y, 
                      X = X,
                      q = 7,
                      seed = 5,  
                      n_steps = 1000, 
                      batch_prop = 1,
                      get_learn_rate = function(i) 1,
                      debug = FALSE, 
                      amortize = FALSE,
                      amortize_in_Y = FALSE,
                      updates = c(Lambda = TRUE, Sigma = TRUE,
                                  Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                  Beta = TRUE, Z = TRUE),
                      get_ELBO_freq = 20,
                      params = result_VI$params)

X_sans_sol <- X[, 1:4]
result_VI_without_Sol <-  get_CAVI(Y = Y, 
                                   X = X_sans_sol,
                                   q = 3,
                                   seed = 5,  
                                   n_steps = 1000, 
                                   batch_prop = 1,
                                   get_learn_rate = function(i) 1,
                                   debug = FALSE, 
                                   amortize = FALSE,
                                   amortize_in_Y = FALSE,
                                   updates = c(Lambda = TRUE, Sigma = TRUE,
                                               Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                               Beta = TRUE, Z = TRUE),
                                   get_ELBO_freq = 20,
                                   params = result_VI_without_Sol$params)

save(result_VI_without_Sol, file = "result_VI_full_CAVI_on_borneo_without_sol.RData")

plot(result_VI_without_Sol$ELBOS)
result_VI <- result_VI_without_Sol
apply(result_VI$params$Lambda$M, 1, var) %>% plot()
round(result_VI$params$Beta$M, 1)
matplot(t(result_VI$params$Beta$M), type = "p", col = 1)


result_VI$params$Beta$M %>% 
  t() %>% 
  matrix(nrow = nrow(.), ncol = ncol(.), 
         dimnames = list(colnames(Y),
                         colnames(X))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Espece") %>% 
  pivot_longer(-c("Espece"),
               names_to = "Variable", values_to = "beta") %>% 
  ggplot(aes(x = beta)) +
  geom_histogram() +
  facet_wrap(~Variable, scales = "free")
