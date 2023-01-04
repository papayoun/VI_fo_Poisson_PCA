rm(list = ls())
library(tidyverse)
library(mixtools)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_PPCA_generating_data_fixed.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_PPCA.txt", sep = ";") %>%
  as.matrix()
X <- read.table("fixed_PPCA.txt", sep = ";") %>%
  as.matrix()
F_x <- ncol(X)

# VI inference functions -----------------------------------------------------

source("utils_PPCA_VI_function.R")

# VI parameters


p <- ncol(Y); q_true <- ncol(Eta_true); n <- nrow(Y)
q <- 7
# on ne prend pas le "vrai" q!
priors <- list(Sigma = list(A = 3, B = 2), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(5, rep(2, q - 1)), 
                           B = 1),
               Beta = list(M = rep(0, F_x),
                           C = rep(0.01, F_x)))
get_result <- function(Y, X, seed, n_steps){
  set.seed(seed)
  if(is.null(X)){
    F_x = 1
  }
  init_params <- list(Lambda = list(M = matrix(rnorm(p * q),
                                               nrow = q, ncol = p),
                                    Cov = array(diag(1, q), 
                                                dim = c(q, q, p))),
                      Beta = list(M = matrix(rnorm(p * F_x),
                                             nrow = F_x, ncol = p),
                                  Cov = array(diag(1, F_x), 
                                              dim = c(F_x, F_x, p))),
                      Eta = list(M = matrix(rnorm(n * q), 
                                            nrow = q, ncol = n),
                                 Cov = array(diag(1, q), 
                                             dim = c(q, q, n))),
                      Sigma = list(A = rep(priors$Sigma$A + n / 2, p),
                                   B = runif(p, 1, 5)),
                      Delta = list(A = priors$Delta$A + 0.5 * p * (q + 1 - (1:q)),
                                   B = runif(q, 1, 10)),
                      Phi = list(A =  matrix(priors$Phi$A + 0.5, p, q),
                                 B = matrix(runif(p * q, 1, 3), p, q))) 
  
  result <- get_CAVI(data_ = Y, X = X, q = q, 
                     n_steps = n_steps, params = init_params,
                     updates = c(Lambda = TRUE, Sigma = TRUE,
                                 Eta = TRUE, Delta = TRUE, 
                                 Phi = TRUE, Beta = !is.null(X)),
                     priors = priors)
  result
}

all_results <- mclapply(1:4, get_result, 
                        Y = Y, X = X,
                        n_steps = 200, 
                        mc.cores = detectCores())


map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 10) %>% 
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line()

params_list <- map(all_results, "params")

best <- params_list[[1]]
t(best$Lambda$M)
Lambda_true

