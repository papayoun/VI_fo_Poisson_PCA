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
init_params <- list(Lambda = list(M = matrix(rnorm(q * p), q, p),
                                  Cov = array(diag(1, q), dim = c(q, q, p))),
                    Eta = list(M = t(Eta_true),
                               Cov = array(diag(1e-4, q), dim = c(q, q, n))),
                    Sigma = list(A = runif(p, 1, 3),
                                 B = runif(p, 1, 3)),
                    Delta = list(A = runif(q, 2, 5) ,
                                 B = rep(1, q)),
                    Phi = list(A = matrix(3/2, p, q),
                               B = matrix(3/2, p, q))) 
  
result <- get_CAVI(data_ = Y, q = q, n_steps = 5, params = NULL,
                   updates = c(Lambda = TRUE, Sigma = TRUE,
                               Eta = FALSE, Delta = TRUE, Phi = TRUE))
params <- result$params
Lambda_true
t(result$params$Lambda$M)
