rm(list = ls())
library(tidyverse)
library(mixtools)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Normal_PPCA_with_covariates.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Normal_PPCA.txt", sep = ";") %>%
  as.matrix()
true_params <- readRDS("experiment_params.rds")
F_x <- ncol(X)
XprimeX = t(X) %*% X

# VI inference functions -----------------------------------------------------

source("utils_Normal_PPCA_VI_functions.R")

# VI parameters


p <- ncol(Y); 
n <- nrow(Y)
q <- ncol(Y) - 1
# on ne prend pas le "vrai" q!
priors <- list(Sigma = list(A = 3, B = 2), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(5, rep(2, q - 1)), 
                           B = 1),
               Beta = list(M = rep(0, F_x),
                           C = rep(0.01, F_x)))

get_result <- function(Y, X, seed, n_steps=n_steps){
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

result=get_result(Y=Y, X=X, seed=1, n_steps=70)

result$ELBOS
all(diff(result$ELBOS[, 2]) > 0)
result$params$Beta$M %>% t() %>% cbind(beta_true)
(round(result$params$Lambda$M[1:4,]*10)/10) %>% as.tibble() %>% t() %>% cbind(Lambda_true)
1:4 %>%
  map(~ get_result(Y=Y, X=X, seed=.x, n_steps=70)) -> 
  all_results

# mclapply ne marche pas????
all_results <- mclapply(1:30,
                        FUN = function(i)
                          get_result(Y = Y, X = X, n_steps = 150, seed = i),
                        mc.cores = detectCores() - 2)

map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 9) %>%
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line() +
  theme(legend.position = "none")

# params_list <- map(all_results, "params")
# 
# best <- params_list[[1]]
# t(best$Lambda$M)
# Lambda_true

(t(result$params$Lambda$M) %*% result$params$Lambda$M + 
  diag(result$params$Sigma$B / result$params$Sigma$A)) %>% 
  cov2cor() %>% 
  round(3)
(true_params$Lambda %*% t(true_params$Lambda) + diag(true_params$sigma2s)) %>% 
  cov2cor() %>% 
  round(3)

result$params$Beta$M
true_params$beta

true_params$Eta %*% t(true_params$Lambda) %>% 
  round(2) %>% 
  head()

t(result$params$Eta$M) %*% result$params$Lambda$M %>% 
  round(2) %>% 
  head()
