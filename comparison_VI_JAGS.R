rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(mvtnorm)
source("utils_generating_data.R") # For true values
source("utils_Poisson_PPCA_VI_functions.R")


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
true_params <- readRDS("experiment_params.rds")

# loading results ---------------------------------------------------------

load("result_VI.RData")
load("result_jags.RData")
head(result_jags)

beta_results_jags <- filter(result_jags,
                            str_detect(Parameter, "beta")) %>%
  mutate(Parameter = droplevels(Parameter))  %>% 
  mutate(method = "JAGS")
table(beta_results_jags$Parameter)


# Tirage selon loi variationnelle -----------------------------------------

beta_results_VI <- map_dfr(1:ncol(Y), 
      function(i){
        mean <- result_VI$params$Beta$M[,i]
        cov <- result_VI$params$Beta$Cov[,,i]
        rmvnorm(1800, mean, cov) %>% 
          matrix(dimnames = list(NULL, paste0("beta[", i, ",", 1:nrow(true_params$beta), "]")),
                 ncol = nrow(true_params$beta)) %>% 
          as.data.frame() %>% 
          rowid_to_column(var = "Iteration") %>% 
          pivot_longer(-Iteration,
                       names_to = "Parameter",
                       values_to = "value")
      }) %>% 
  mutate(method = "VI")

truebetas = data.frame(Parameter=map(1:ncol(true_params$beta), 
                                     function(i) {paste0("beta[", i, ",", 1:nrow(true_params$beta), "]")}) %>% 
                         unlist(), intercept= as.vector(true_params$beta))

bind_rows(beta_results_VI,
          beta_results_jags) %>% 
  ggplot() + 
  aes(x = value) +
  geom_density(aes(color = method)) +
  geom_vline(data=truebetas, mapping = aes(xintercept = intercept)) +
  geom_point(data=truebetas,mapping= aes(x=intercept, y=0))+
  facet_wrap(~Parameter, ncol = 2)


bind_rows(beta_results_VI,
          beta_results_jags) %>% 
  filter(str_detect(Parameter, "beta\\[2")) %>% 
  pivot_wider(values_from = "value", names_from = "Parameter") %>% 
  ggplot(aes(x = `beta[2,1]`, y = `beta[2,2]`)) +
  geom_point(aes(color = method))+
  geom_point(data=truebetas %>% 
               filter(str_detect(Parameter, "beta\\[2")) %>% 
               pivot_wider(values_from = "intercept", names_from = "Parameter"), color="black",shape= 18, size=10)

# loading results ---------------------------------------------------------


sigma_results_jags <- filter(result_jags,
                            str_detect(Parameter, "preciE")) %>%
  mutate(Parameter = droplevels(Parameter))  %>% 
  mutate(method = "JAGS")



sigma_results_VI <- map_dfr(1:ncol(Y), 
                           function(i){
                             A <- result_VI$params$Sigma$A[i]
                             B <- result_VI$params$Sigma$B[i]
                             rgamma(1800, A, B) %>% 
                               matrix(dimnames = list(NULL, paste0("preciE[", i, "]")),
                                      ncol = 1) %>% 
                               as.data.frame() %>% 
                               rowid_to_column(var = "Iteration") %>% 
                               pivot_longer(-Iteration,
                                            names_to = "Parameter",
                                            values_to = "value")
                           }) %>% 
  mutate(method = "VI")

bind_rows(sigma_results_VI ,
          sigma_results_jags) %>% 
  group_by(Parameter, method) %>% 
  summarise(msigmamoins2= mean(value))
1/true_params$sigma2s 

bind_rows(sigma_results_VI,
          sigma_results_jags) %>% 
  ggplot() + 
  aes(x = value) +
  geom_density(aes(color = Parameter)) +
  facet_wrap(~method, ncol = 2, scales = "free")


dim(Y)[2]

# Lambdas -----------------------------------------------------------------

n_samples <- 200
Lambdas_jags <- filter(result_jags, 
       str_detect(Parameter, "lambda")) %>% 
  group_by(Chain, Iteration) %>% 
  group_map(function(x, g){
    pull(x, value) %>% 
      matrix(nrow = ncol(Y), byrow = TRUE)
  }) %>% 
  {.[1:n_samples]}
Etas_jags <- filter(result_jags, 
                       str_detect(Parameter, "beta", negate = TRUE),
                    str_detect(Parameter, "eta")) %>% 
  group_by(Chain, Iteration) %>% 
  group_map(function(x, g){
    pull(x, value) %>% 
      matrix(nrow = nrow(Y), byrow = TRUE)
  })%>% 
  {.[1:n_samples]}

Betas_jags <- filter(result_jags, 
                    str_detect(Parameter, "beta")) %>% 
  group_by(Chain, Iteration) %>% 
  group_map(function(x, g){
    pull(x, value) %>% 
      matrix(nrow = ncol(Y), byrow = TRUE)
  })%>% 
  {.[1:n_samples]}
Lambda_Eta_jags <- pmap(.l= list(x = Etas_jags, y = Lambdas_jags, z = Betas_jags),
                        .f = function(x, y, z){
  x %*% t(y) + X %*% t(z)
})

Etas_VI <- rerun(length(Lambda_Eta_jags),
                 sapply(1:nrow(Y),
       function(i){
         mixtools::rmvnorm(n = 1, 
                 mu = result_VI$params$Eta$M[,i], 
                 sigma = result_VI$params$Eta$Cov[,,i])
       }) %>% t())
Lambdas_VI <- rerun(length(Lambda_Eta_jags),
                 sapply(1:ncol(Y),
                        function(j){
                          mixtools::rmvnorm(n = 1, 
                                           mu = result_VI$params$Lambda$M[,j], 
                                           sigma = result_VI$params$Lambda$Cov[,,j])
                        }) %>%
                   t())
Betas_VI <- rerun(length(Lambda_Eta_jags),
                 sapply(1:ncol(Y),
                        function(i){
                          mixtools::rmvnorm(n = 1, 
                                            mu = result_VI$params$Beta$M[,i], 
                                            sigma = result_VI$params$Beta$Cov[,,i])
                        }) %>% t())
Lambda_Eta_VI <- pmap(.l= list(x = Etas_VI, y = Lambdas_VI, z = Betas_VI),
                      .f = function(x, y, z){
                        x %*% t(y) + X %*% t(z)
                      })

Lambda_Eta_VI_df <- do.call(what = function(...) abind(..., along = 3), 
                            args = Lambda_Eta_VI) %>% 
  format_array(param_name = "Lambda_Eta") %>% 
  mutate(method = "VI")
Lambda_Eta_jags_df <- do.call(what = function(...) abind(..., along = 3), 
                            args = Lambda_Eta_jags) %>% 
  format_array(param_name = "Lambda_Eta") %>% 
  mutate(method = "jags")
format_array(array(Y, dim = c(nrow(Y), ncol(Y), 1)),
              param_name = "Y")
bind_rows(Lambda_Eta_VI_df,
          Lambda_Eta_jags_df) %>% 
  separate(Parameter, into = c("Nom", "Indice"), sep = "Lambda_Eta") %>% 
  ggplot(aes(x = Estimate, color = method)) + 
  geom_density() +
  facet_wrap(Indice~., nrow = 30) +
  theme(strip.background = element_blank(), 
        strip.text = element_blank()) +
  geom_vline(data = format_array(array(Y, dim = c(nrow(Y), ncol(Y), 1)),
                                 param_name = "Y") %>% 
               separate(Parameter, into = c("Nom", "Indice"), sep = "Y"),
             aes(xintercept = log(Estimate + .01)))
site <- sample(1:nrow(Y), size = 1)
species <- sample(1:ncol(Y), size = 1)

tibble(value = map_dbl(Lambda_Eta_VI, function(x) x[site, species]),
       method = "VI") %>% 
  bind_rows(tibble(value = map_dbl(Lambda_Eta_jags, function(x) x[site, species]),
                   method = "Jags")) %>% 
  ggplot(aes(x = value, color = method)) +
  geom_density() +
  geom_vline(xintercept = log(Y[site, species] + .1)) +
  labs(title = site, subtitle = species)



Lambdas_VI
apply(result_VI$params$Lambda$M, 1, var)
result_VI$params$Phi$A / result_VI$params$Phi$B
plot(1 / cumprod(result_VI$params$Delta$A / result_VI$params$Delta$B))
lambda_results_jags <- filter(result_jags,
                            str_detect(Parameter, "lambda")) %>%
  mutate(Parameter = droplevels(Parameter))  %>% 
  mutate(method = "JAGS")
lambda_results_jags %>% 
  group_by(Parameter) %>% 
  summarise(Lambda_jags_mean = mean(value)) %>% pull() %>% matrix(nc=7,byrow = T)->Lambda_jags_mean


Lambdas_jags %>% 
  sapply(function(x) apply(x, 2, var)) %>% 
  t() %>% 
  boxplot()
Lambdas_VI %>% 
  sapply(function(x) apply(x, 2, var)) %>% 
  t() %>% 
  boxplot()

# Eta Lambdas -------------------------------------------------------------

result_jags$

