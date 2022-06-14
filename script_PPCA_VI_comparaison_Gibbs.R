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
result_VI$Eta<-array(0,dim = c(n,q,n_samples))
for(i in 1:n){result_VI$Lambda[i,1:q,1:n_samples]=
  t(mvtnorm::rmvnorm(n_samples,result_VI_coef$Lambda$M[1:q,i],
                     result_VI_coef$Lambda$Cov[1:q,1:q,i]))}
result_VI$Sigma<-matrix(0,nr=p,nc=n_samples)
for(s in 1:n_samples){result_VI$Sigma[1:p,s]=1/rgamma(p,result_VI_coef$Sigma$A[1:p],
                             result_VI_coef$Sigma$B[1:p])}
result_VI$deltas<-matrix(0,nr=q,nc=n_samples)
for(s in 1:n_samples){result_VI$deltas[1:q,s]=rgamma(q,
                                result_VI_coef$Delta$A[1:q],
                                result_VI_coef$Delta$B[1:q])}

###########

#prep result_VI
LL_t_VI <- mclapply(1:n_samples, function(i){
  result_VI$Lambda[,,i] %*% t(result_VI$Lambda[,,i]) 
}) %>% 
  abind(along = 3)
LL_t_hat_VI_df <- format_array(LL_t_VI, "LL_t") 
LL_t_hat_VI_df<-LL_t_hat_VI_df %>% mutate(method="VI")

#prep result_gibbs
LL_t_gibbs <- mclapply(1:n_samples, function(i){
  result_gibbs$Lambda[,,i] %*% t(result_gibbs$Lambda[,,i]) 
}) %>% 
  abind(along = 3)
LL_t_hat_gibbs_df <- format_array(LL_t_gibbs, "LL_t") 
LL_t_hat_gibbs_df <-LL_t_hat_gibbs_df %>% mutate(method="Gibbs")
#verité et les deux methode
LL_t_true <- Lambda_true %*% t(Lambda_true)
LL_t_hat_df <-LL_t_hat_gibbs_df %>% bind_rows(LL_t_hat_VI_df)
# une densité par parametre
  ggplot(LL_t_hat_df) + 
    aes(x = Estimate, color=method) +
    geom_density() +
    facet_wrap(~Parameter, nrow =  nrow(Lambda_true), labeller = label_parsed,
               scales = "free_y") +
    geom_vline(data = tibble(Parameter = unique(LL_t_hat_df$Parameter),
                             Truth = as.numeric(t(LL_t_true))),
               aes(xintercept = Truth),
               color = "black")
  
  Lambda_hat_gibbs_df <- format_array(result_gibbs$Lambda,
                                "Lambda") %>% mutate(method="Gibbs")
  Lambda_hat_VI_df <- format_array(result_VI$Lambda,
                                      "Lambda") %>% mutate(method="VI")
  Lambda_hat_df<-Lambda_hat_gibbs_df %>% bind_rows(Lambda_hat_VI_df)
  Lambda_hat_df %>% 
    filter((str_detect(Parameter, "-1]") | str_detect(Parameter, "-2]"))) %>% 
    separate(Parameter, into = c("species", "feature"), sep = "-") %>% 
    mutate(feature = paste0("feature_", str_remove(feature , "]")),
           species = str_replace(species, "Lambda\\[", "species_")) %>%
    pivot_wider(names_from = "feature", values_from = "Estimate") %>% 
    ggplot(aes(x = feature_1, y = feature_2, color = species, size=method, shape=method)) +
    geom_point() 
  ########################Comparaison des sigma^2
  Sigma_gibbs_df<- format_matrix(matrix_=result_gibbs$Sigma,
                                 param_name="s",
                                 suffix_ = "2") %>% 
    mutate(method="Gibbs")
  Sigma_VI_df<- format_matrix(matrix_=result_VI$Sigma,
                              param_name="s",
                              suffix_ = "2")%>% 
    mutate(method="VI")
  Sigma_df<-Sigma_VI_df %>% bind_rows(Sigma_gibbs_df)
  Sigma_df %>% ggplot(aes(x=Parameter,y=Estimate,color=method))+
    geom_boxplot()
  ########################Comparaison des delta##########
  Delta_gibbs_df<- format_matrix(matrix_=result_gibbs$deltas,
                                 param_name="delta") %>% 
    mutate(method="Gibbs")
  Delta_VI_df<- format_matrix(matrix_=result_VI$deltas,
                              param_name="delta")%>% 
    mutate(method="VI")
  Delta_df<-Delta_VI_df %>% bind_rows(Delta_gibbs_df)
  Delta_df %>% ggplot(aes(x=Parameter,y=Estimate,color=method))+
    geom_boxplot()
  