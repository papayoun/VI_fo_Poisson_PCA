
rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(pacman) # For progress bar
library(torch)
library(fastDummies)

# Data --------------------------------------------------------------------

source("utils_load_main_application_data.R")

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_functions.R")

# Inference with complete covariates --------------------------------------

# # X_complete <- select(X_full, -Site) %>% 
# #   as.matrix()
# 
# # result_VI <- NULL
# # if (file.exists("result_VI_full_CAVI_on_borneo.RData"))
# #   load("result_VI_full_CAVI_on_borneo.RData")
# # result_VI <- get_CAVI(Y = Y_complete, 
# #                       X = X_complete,
# #                       q = 7,
# #                       seed = 5,  
# #                       n_steps = 1000, 
# #                       batch_prop = 1,
# #                       get_learn_rate = function(i) 1,
# #                       debug = FALSE, 
# #                       amortize = FALSE,
# #                       amortize_in_Y = FALSE,
# #                       updates = c(Lambda = TRUE, Sigma = TRUE,
# #                                   Eta = TRUE, Delta = TRUE, Phi = TRUE, 
# #                                   Beta = TRUE, Z = TRUE),
# #                       get_ELBO_freq = 50,
# #                       params = result_VI$params)
# # save(result_VI, file = "result_VI_full_CAVI_on_borneo.RData")
# # 
# # rm(list = ls(pattern = "result_VI"))
# # rm(list = ls(pattern = "X_complete"))
# 

# Without Sol -------------------------------------------------------------


result_VI_without_sol <- NULL
if (file.exists("result_VI_full_CAVI_on_borneo_without_sol.RData"))
  load("result_VI_full_CAVI_on_borneo_without_sol.RData")

X_without_sol <- X_full %>% 
  select_at(vars(-starts_with("Sol"), -Site)) %>% 
  as.matrix() %>% 
  cbind(Intercept = 1, .)
result_VI_without_sol <-  get_CAVI(Y = Y_complete, 
                                   X = X_without_sol,
                                   q = 7,
                                   seed = 5,  
                                   n_steps = 200, 
                                   batch_prop = 1,
                                   get_learn_rate = function(i) 1,
                                   debug = FALSE, 
                                   amortize = FALSE,
                                   amortize_in_Y = FALSE,
                                   updates = c(Lambda = TRUE, Sigma = TRUE,
                                               Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                               Beta = TRUE, Z = TRUE),
                                   get_ELBO_freq = 40,
                                   params = result_VI_without_sol$params)


save(result_VI_without_sol, file = "result_VI_full_CAVI_on_borneo_without_sol.RData")

rm(list = ls(pattern = "result_VI_without_sol"))
rm(list = ls(pattern = "X_complete"))
