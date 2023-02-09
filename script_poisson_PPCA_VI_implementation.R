rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_Poisson_PPCA_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
X_data <- read.table("fixed_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_function.R")

# VI parameters


p <- ncol(Y); q_true <- ncol(Eta_true); n <- nrow(Y)
q <- 7
# on ne prend pas le "vrai" q!
priors <- list(Sigma = list(A = 3, B = 2), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(5, rep(2, q - 1)), 
                           B = 1),
               Beta = list(M = rep(0, ncol(X_data)),
                           C = rep(0.01, ncol(X_data))))

get_result <- function(Y, X, seed, n_steps=n_steps, 
                       updates = c(Lambda = TRUE, Sigma = TRUE,
                                   Eta = TRUE, Delta = TRUE, 
                                   Phi = TRUE, Beta = !is.null(X), Z = TRUE), 
                       debug = FALSE){
  set.seed(seed)
  if(is.null(X)){
    F_x = 1
  }
  else{
    F_x <- ncol(X)
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
                                 B = matrix(runif(p * q, 1, 3), p, q)),
                      Z = list(M = log(Y + 1),
                               S2 = matrix(.1, nrow = nrow(Y),
                                           ncol = ncol(Y)))) 
  result <- get_CAVI(data_ = Y, X = X, q = q, 
                     n_steps = n_steps, params = init_params,
                     updates = updates,
                     priors = priors, debug = debug)
  result
}

result=get_result(Y=Y, 
                  X = X_data, 
                  seed = 1, 
                  n_steps = 50, 
                  updates = c(Lambda = TRUE, Sigma = TRUE,
                              Eta = TRUE, Delta = TRUE, 
                              Phi = TRUE, Beta = TRUE, Z = TRUE),
                  debug = TRUE)
plot(result$ELBOS)
t(result$params$Beta$M)
beta_true
all(diff(result$ELBOS[, 2]) > 0)
result$params$Beta$M %>% t() %>% cbind(beta_true)

1:4 %>%
  map(~ get_result(Y = Y, X = X_data, seed = .x, n_steps = 70)) -> 
  all_results

all_results <- mclapply(1:30,
                        FUN = function(i)
                          get_result(Y = Y, X = X_data, n_steps = 150, seed = i),
                        mc.cores = detectCores() - 2)

map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 120) %>%
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line() +
  theme(legend.position = "none")

params_list <- map(all_results, "params")

best <- params_list[[1]]
t(best$Beta$M)
beta_true
# t(best$Lambda$M)
# Lambda_true

