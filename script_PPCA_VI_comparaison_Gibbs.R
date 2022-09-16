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

q_guess <- 7
# on ne prend pas le "vrai" q!
priors <- list(Sigma = list(A = 3, B = 2), 
               Phi = list(A = 3/2, B = 3/2),
               Delta= list(A = c(5, rep(2, q_guess - 1)), 
                           B = 1))
get_result <- function(seed, n_steps){
  set.seed(seed)
  init_params <- list(Lambda = list(M = matrix(rnorm(p * q_guess),
                                               nrow = q_guess, ncol = p),
                                    Cov = array(diag(1, q_guess), 
                                                dim = c(q_guess, q_guess, p))),
                      Eta = list(M = matrix(rnorm(n * q_guess), 
                                            nrow = q_guess, ncol = n),
                                 Cov = array(diag(1, q_guess), 
                                             dim = c(q_guess, q_guess, n))),
                      Sigma = list(A = rep(priors$Sigma$A + n / 2, p),
                                   B = runif(p, 1, 5)),
                      Delta = list(A = priors$Delta$A + 0.5 * p * (q_guess + 1 - (1:q_guess)),
                                   B = runif(q_guess, 1, 10)),
                      Phi = list(A =  matrix(priors$Phi$A + 0.5, p, q_guess),
                                 B = matrix(runif(p * q_guess, 1, 3), p, q_guess))) 
  
  result <- get_CAVI(data_ = Y, q = q_guess, 
                     n_steps = n_steps, params = init_params,
                     updates = c(Lambda = TRUE, Sigma = TRUE,
                                 Eta = TRUE, Delta = TRUE, 
                                 Phi = TRUE),
                     priors = priors)
  result
}

all_results <- mclapply(1:16, get_result, 
                        n_steps = 2000, 
                        mc.cores = detectCores())
best_index <- map_dbl(all_results, 
                      function(res) 
                        res$ELBOS[nrow(res$ELBOS), "ELBO"]) %>% 
  which.max()

map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration > 100) %>% 
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line()

params_list <- map(all_results, "params")

result_VI_coef<- params_list[[best_index]]

####################
n_steps <- 10000
n_burn <- 50
n_thin <- 50
k_tilde <- q_guess

result_gibbs <- get_gibbs_sample(data_ = Y, n_steps = n_steps, k_tilde = k_tilde, 
                                 burn = n_burn, thin = n_thin)

result_gibbs$tausPhis <- lapply(1:p,
                                function(j)
                                  result_gibbs$Phi[j,,] * result_gibbs$taus) %>% 
  do.call(what = function(...) abind(..., along = 3)) %>% 
  aperm(perm = c(3, 1, 2))
#####################
n_samples <-n_steps/n_thin
result_VI<-list()
# Faisons une structure similaire à Gibbs
result_VI$Lambda<-array(0,dim = c(p,q_guess,n_samples))
for(j in 1:p){result_VI$Lambda[j,1:q_guess,1:n_samples]=
  t(mvtnorm::rmvnorm(n_samples,result_VI_coef$Lambda$M[1:q_guess,j],
                     result_VI_coef$Lambda$Cov[1:q_guess,1:q_guess,j]))}
result_VI$Eta<-array(0,dim = c(n,q_guess,n_samples))
for(i in 1:n){result_VI$Eta[i,1:q_guess,1:n_samples]=
  t(mvtnorm::rmvnorm(n_samples,result_VI_coef$Eta$M[,i],
                     result_VI_coef$Eta$Cov[,,i]))}
result_VI$Sigma<-matrix(0,nr=p,nc=n_samples)
for(s in 1:n_samples){result_VI$Sigma[1:p,s]=1/rgamma(p,result_VI_coef$Sigma$A[1:p],
                                                      result_VI_coef$Sigma$B[1:p])}
result_VI$deltas<-matrix(0,nr=q_guess,nc=n_samples)
for(s in 1:n_samples){result_VI$deltas[1:q_guess,s]=rgamma(q_guess,
                                                     result_VI_coef$Delta$A[1:q_guess],
                                                     result_VI_coef$Delta$B[1:q_guess])}
result_VI$taus <- apply(result_VI$deltas, 2, cumprod)
result_VI$Phi <- rerun(n_samples,
                       sapply(1:q_guess, 
                              function(h) 
                                rgamma(p, 
                                       result_VI_coef$Phi$A[,h],
                                       result_VI_coef$Phi$B[,h]))) %>% 
  do.call(what = function(...) abind(..., along = 3))
result_VI$tausPhis <- lapply(1:p,
                             function(j)
                               result_VI$Phi[j,,] * result_VI$taus) %>% 
  do.call(what = function(...) abind(..., along = 3)) %>% 
  aperm(perm = c(3, 1, 2))
###########
EtaEta_t_VI <- mclapply(1:n_samples, function(i){
  t(result_VI$Eta[,,i]) %*% (result_VI$Eta[,,i]) 
}) %>% 
  abind(along = 3)

EtaEta_t_gibbs <- mclapply(1:n_samples, function(i){
  t(result_gibbs$Eta[,,i]) %*% (result_gibbs$Eta[,,i]) 
}) %>% 
  abind(along = 3)

EtaEta_t_df <- format_array(EtaEta_t_VI, "EE_t") %>%
  mutate(method="VI") %>% 
  bind_rows(format_array(EtaEta_t_gibbs, "EE_t") %>%
              mutate(method="Gibbs"))
ggplot(EtaEta_t_df) + 
  aes(x = Estimate, color=method) +
  geom_density() +
  facet_wrap(~Parameter, 
             ncol =  q_guess, 
             labeller = label_parsed,
             scales = "free_y")

result_VI_coef$Eta$Cov[,,37]


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
  ggplot(aes(x = feature_1, y = feature_2, 
             color = species)) +
  facet_wrap(~method) +
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


tausPhis_df <- format_array(result_gibbs$tausPhis,
                            param_name = "tausPhis") %>% 
  mutate(method = "Gibbs") %>% 
  bind_rows(format_array(result_VI$tausPhis,
                         param_name = "tausPhis") %>% 
              mutate(method = "VI"))

tausPhis_df %>% 
  filter(str_detect(Parameter, "\\[5-")) %>% 
  ggplot(aes(x = Parameter, 
             y = Estimate, fill = method)) +
  geom_boxplot() +
  facet_wrap(~Parameter, scales = "free") +
  scale_y_continuous(trans = "log10")

Lambda_hat_df %>% 
  filter(str_detect(Parameter, "\\[5-")) %>% 
  ggplot(aes(x = Parameter, 
             y = Estimate, fill = method)) +
  geom_boxplot() +
  facet_wrap(~Parameter, scales = "free")

Lambda_hat_df %>% 
  separate(Parameter, 
           into = c("ligne", "colonne"),
           sep = "-") %>% 
  group_by(colonne, iteration, method) %>% 
  summarise(variance = var(Estimate)) %>% 
  ggplot(aes(x = colonne, y = variance, fill = method)) +
  geom_boxplot()
