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
priors <- list(Sigma = list(A = 3, B = 2), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(5, rep(2, q - 1)), 
                           B = 1))
init_params <- list(Lambda = list(M = t(Lambda_true),
                                  Cov = array(diag(1e-4, q), 
                                              dim = c(q, q, p))),
                    Eta = list(M = matrix(rnorm(n * q), 
                                          nrow = q, ncol = n),
                               Cov = array(diag(1, q), 
                                           dim = c(q, q, n))),
                    Sigma = list(A = rep(2, p),
                                 B = rep(3, p)),
                    Delta = list(A = sort(runif(q, 2, 5), decreasing = TRUE ,
                                 B = rep(1, q)),
                    Phi = list(A = matrix(3/2, p, q),
                               B = matrix(3/2, p, q))) 
  
result <- get_CAVI(data_ = Y, q = q, 
                   n_steps = 10, params = init_params,
                   updates = c(Lambda = TRUE, Sigma = TRUE,
                               Eta = TRUE, Delta = TRUE, 
                               Phi = TRUE),
                   priors = priors)
Lambda_est <- t(result$params$Lambda$M) 
Lambda_true %*% t(Lambda_true)
round(Lambda_est %*% t(Lambda_est), 3)
params <- result$params
head(Eta_true)
head(t(result$params$Eta$M))
result$ELBOS %>% 
  ggplot() + 
  aes(x = iteration, y = ELBO) + 
  geom_line()
any(diff(result$ELBOS$ELBO) < 0)
