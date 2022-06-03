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
get_result <- function(seed, n_steps){
  set.seed(seed)
  init_params <- list(Lambda = list(M = matrix(rnorm(p * q),
                                               nrow = q, ncol = p),
                                    Cov = array(diag(1, q), 
                                                dim = c(q, q, p))),
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
  
  result <- get_CAVI(data_ = Y, q = q, 
                     n_steps = n_steps, params = init_params,
                     updates = c(Lambda = FALSE, Sigma = TRUE,
                                 Eta = TRUE, Delta = TRUE, 
                                 Phi = TRUE),
                     priors = priors)
  result
}

all_results <- mclapply(1:8, get_result, n_steps = 50, mc.cores = 4)
map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 9) %>% 
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line()
map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  group_by(Replicate) %>% 
  filter(iteration > 9)

params <- result$params
head(Eta_true)
head(t(result$params$Eta$M))
result$ELBOS %>% 
  ggplot() + 
  aes(x = iteration, y = ELBO) + 
  geom_line()
any(diff(result$ELBOS$ELBO) < 0)
