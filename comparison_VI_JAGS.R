rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_generating_data.R") # For true values
source("utils_Poisson_PPCA_VI_functions.R") 


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()


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

beta_results_VI <- map_dfr(1:10, 
      function(i){
        mean <- result_VI$params$Beta$M[,i]
        cov <- result_VI$params$Beta$Cov[,,i]
        rmvnorm(1800, mean, cov) %>% 
          matrix(dimnames = list(NULL, paste0("beta[", i, ",", 1:3, "]")),
                 ncol = 3) %>% 
          as.data.frame() %>% 
          rowid_to_column(var = "Iteration") %>% 
          pivot_longer(-Iteration,
                       names_to = "Parameter",
                       values_to = "value")
      }) %>% 
  mutate(method = "VI")

truebetas= data_frame(Parameter=map(1:10, function(i) {paste0("beta[", i, ",", 1:3, "]")}) %>% unlist(), intercept= as.vector(beta_true))

bind_rows(beta_results_VI,
          beta_results_jags) %>% 
  ggplot() + 
  aes(x = value) +
  geom_density(aes(color = method)) +
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



sigma_results_VI <- map_dfr(1:10, 
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
          sigma_results_jags) %>% group_by(Parameter, method) %>% summarise(msigmamoins2= mean(value))
sigma_m2s_true 

bind_rows(sigma_results_VI,
          sigma_results_jags) %>% 
  ggplot() + 
  aes(x = value) +
  geom_density(aes(color = method)) +
  facet_wrap(~Parameter, ncol = 1)




# Lambdas -----------------------------------------------------------------


lambda_results_jags <- filter(result_jags,
                            str_detect(Parameter, "lambda")) %>%
  mutate(Parameter = droplevels(Parameter))  %>% 
  mutate(method = "JAGS")
lambda_results_jags %>% 
  group_by(Parameter) %>% 
  summarise(Lambda_jags_mean = mean(value)) %>% pull() %>% matrix(nc=7,byrow = T)->Lambda_jags_mean
Lambda_true
Lambda_jags_mean 
abs(Lambda_jags_mean)>0.25
(abs(Lambda_jags_mean)>0.25)*round(Lambda_jags_mean,3)
Lambda_true

