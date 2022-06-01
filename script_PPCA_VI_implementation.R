rm(list = ls())
library(tidyverse)
library(mixtools)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_PPCA_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_PPCA.txt", sep = ";") %>%
  as.matrix()

# VI inference functions -----------------------------------------------------

source("utils_PPCA_VI_function.R")

# VI parameters


p <- ncol(Y); q <- ncol(Eta_true); n <- nrow(Y)
priors <- list(Sigma = list(A = 1, B = 3), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(2, rep(3, q - 1)), 
                           B = 1))
init_params <- list(Lambda = list(M = matrix(rnorm(q * p), q, p),
                                  Cov = array(diag(1, q), dim = c(q, q, p))),
                    Eta = list(M = t(Eta_true),
                               Cov = array(diag(1e-4, q), dim = c(q, q, n))),
                    Sigma = list(A = sigma_m2s_true^2 * 1e3,
                                 B = sigma_m2s_true * 1e3),
                    Delta = list(A = runif(q, 2, 5) ,
                                 B = rep(1, q)),
                    Phi = list(A = matrix(3/2, p, q),
                               B = matrix(3/2, p, q))) 
  
result <- get_CAVI(data_ = Y, q = q, n_steps = 10, params = init_params,
                   updates = c(Lambda = TRUE, Sigma = FALSE,
                               Eta = TRUE, Delta = TRUE, Phi = TRUE),
                   priors = priors)
params <- result$params
Lambda_true
t(result$params$Lambda$M)
result$ELBOS %>% 
  ggplot() + 
  aes(x = iteration, y = ELBO) + 
  geom_line()
