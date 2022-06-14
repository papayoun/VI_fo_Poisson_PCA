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
source("utils_PPCA_gibbs_function.R")

# VI parameters


p <- ncol(Y); q_true <- ncol(Eta_true); n <- nrow(Y)
q<-7
# on ne prend pas le "vrai" q!
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
                     updates = c(Lambda = TRUE, Sigma = TRUE,
                                 Eta = TRUE, Delta = TRUE, 
                                 Phi = TRUE),
                     priors = priors)
  result
}

all_results <- mclapply(1:4, get_result, 
                        n_steps = 2000, 
                        mc.cores = detectCores())


map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 10) %>% 
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line()

params_list <- map(all_results, "params")

best <- params_list[[1]]
t(best$Lambda$M)
Lambda_true
result_VI_coef<- params_list[[1]]

####################
n_steps <- 10000
n_burn <- 50
n_thin <- 50
k_tilde <- q

result_gibbs <- get_gibbs_sample(data_ = Y, n_steps = n_steps, k_tilde = k_tilde, 
                           burn = n_burn, thin = n_thin)
#####################
n_samples <-n_steps/n_thin
result_VI<-list()
# Faisons une structure similaire à Gibbs
result_VI$Lambda<-array(0,dim = c(p,q,n_samples))
for(j in 1:p){result_VI$Lambda[j,1:q,1:n_samples]=
t(mvtnorm::rmvnorm(n_samples,result_VI_coef$Lambda$M[1:q,j],
                   result_VI_coef$Lambda$Cov[1:q,1:q,j]))}
###########
result=result_VI
LL_t <- mclapply(1:n_samples, function(i){
  result$Lambda[,,i] %*% t(result$Lambda[,,i]) 
}) %>% 
  abind(along = 3)

LL_t_hat_df <- format_array(LL_t, "LL_t") 
LL_t_true <- Lambda_true %*% t(Lambda_true)

LL_t_hat_df %>% 
  ggplot(aes(x = iteration, y = Estimate, color = Parameter)) +
  geom_line()
  ggplot(LL_t_hat_df) + 
    aes(x = Estimate) +
    geom_density() +
    facet_wrap(~Parameter, nrow =  nrow(Lambda_true), labeller = label_parsed,
               scales = "free_y") +
    geom_vline(data = tibble(Parameter = unique(LL_t_hat_df$Parameter),
                             Truth = as.numeric(t(LL_t_true))),
               aes(xintercept = Truth),
               color = "red")
