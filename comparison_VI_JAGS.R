rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(mvtnorm)
# source("utils_generating_data.R") # For true values
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
beta_results_jags %>% 
  group_by(Parameter) %>% 
  summarise(est = mean(value)) %>% 
  pull(est) %>% 
  matrix(nrow = 2) %>% 
  round(2)
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


# result_jags %>% dplyr::filter(str_detect(Parameter,"beta")) %>% 
#   group_by(Parameter) %>% summarize(value=mean(value)) ->tt

# i=1
# j=1
# paste0("beta[",i,",",j,"]")
# tt %>% dplyr::filter(str_detect(Parameter,"beta[1,1]")) %>% pull(value)
# 
# nbi=dim(result_VI$params$Beta$M)[1]
# nbj=dim(result_VI$params$Beta$M)[2]
# bb<-result_VI$params$Beta$M*0
# for (j in 1:nbj) {
#   for (i in 1:nbi) {bb[i,j]<-tt}}

# Structure pour approximer la conjointe
params_Jags = list()
# Z
nbi=dim(result_VI$params$Z$M)[1]
nbj=dim(result_VI$params$Z$M)[2]
tabZ <- result_jags %>% dplyr::filter(str_detect(Parameter,"Z")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain)) 

# Moyenne et variance des Z -----------------------------------------------------------

tabZ %>% apply(2,mean) %>% matrix(nrow = nbi, ncol = nbj, byrow=T) -> params_Jags$Z$M
tabZ %>% apply(2,var) %>% matrix(nrow = nbi, ncol = nbj, byrow=T) -> params_Jags$Z$S2
rm(tabZ)

# beta

nbi=dim(result_VI$params$Beta$M)[1]
nbj=dim(result_VI$params$Beta$M)[2]
result_jags %>% dplyr::filter(str_detect(Parameter,"beta")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabBeta 
tabBeta %>% apply(2,mean) %>% 
  matrix(nrow = nbj,ncol = nbi,byrow=T) %>% t(.) -> params_Jags$Beta$M
params_Jags$Beta$Cov<- 0*result_VI$params$Beta$Cov
for(j in 1:nbj){
  pattern<-paste0("[",j,",")
  tabBeta %>% select(contains(pattern)) %>% var -> params_Jags$Beta$Cov[,,j]
}
rm(tabBeta,pattern)
# phi
nbi=dim(result_VI$params$Phi$A)[1]
nbj=dim(result_VI$params$Phi$A)[2]
result_jags %>% dplyr::filter(str_detect(Parameter,"phi")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabPhi
tabPhi %>% apply(2,mean) %>% 
  matrix(nrow = nbj,ncol = nbi,byrow=F) %>% t(.)->MM
tabPhi %>% apply(2,var) %>% 
  matrix(nrow = nbj,ncol = nbi,byrow=F) %>% t(.)->VV
params_Jags$Phi$A <- MM^2/VV
params_Jags$Phi$B <- MM/VV
rm(tabPhi,MM,VV)
# delta
nbi=dim(result_VI$params$Delta$A)
result_jags %>% dplyr::filter(str_detect(Parameter,"delta")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabDelta
tabDelta %>% apply(2,mean) ->MM
tabDelta %>% apply(2,var) ->VV
params_Jags$Delta$A<- MM^2/VV
params_Jags$Delta$B<- MM/VV
rm(tabDelta,MM,VV)
# sigma
nbi=dim(result_VI$params$Sigma$A)
result_jags %>% dplyr::filter(str_detect(Parameter,"preciE")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabSigma
tabSigma %>% apply(2,mean) ->MM
tabSigma %>% apply(2,var) ->VV
params_Jags$Sigma$A<- MM^2/VV
params_Jags$Sigma$B<- MM/VV
rm(tabSigma,MM,VV)
#Eta
nbi=dim(result_VI$params$Eta$M)[1]
nbj=dim(result_VI$params$Eta$M)[2]
result_jags %>% dplyr::filter(str_detect(Parameter,"beta",negate=T)) %>% 
  dplyr::filter(str_detect(Parameter,"eta")) %>%   pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabEta 
tabEta %>% apply(2,mean) %>% 
  matrix(nrow = nbi,ncol = nbj,byrow=F) -> params_Jags$Eta$M
params_Jags$Eta$Cov<- 0*result_VI$params$Eta$Cov
for(j in 1:nbj){
  pattern<-paste0("[",j,",")
  tabEta %>% select(contains(pattern)) %>% var -> params_Jags$Eta$Cov[,,j]
}
rm(tabEta,pattern)
# lambda
nbi=dim(result_VI$params$Lambda$M)[1]
nbj=dim(result_VI$params$Lambda$M)[2]
result_jags %>% dplyr::filter(str_detect(Parameter,"lambda")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain))  ->tabLambda 
tabLambda %>% apply(2,mean) %>% 
  matrix(nrow = nbj,ncol = nbi,byrow=T) %>% t(.) -> params_Jags$Lambda$M
params_Jags$Lambda$Cov<- 0*result_VI$params$Lambda$Cov
for(j in 1:nbj){
  pattern<-paste0("[",j,",")
  tabLambda %>% select(contains(pattern)) %>% var -> params_Jags$Lambda$Cov[,,j]
}
rm(tabLambda,pattern)

n=result_VI$params$n
p=result_VI$params$p
q=result_VI$params$q
F_x=result_VI$params$F_x
params_Jags$n=result_VI$params$n
params_Jags$p=result_VI$params$p
params_Jags$q=result_VI$params$q
params_Jags$F_x=result_VI$params$F_x
### CALCUL ELBO

priors = list(Sigma = list(A = 3, B = 2), 
              Phi = list(A = 3/2, B = 3/2),
              Delta= list(A = rep(3, q), 
                          B = 1))
priors$Beta = list(M = matrix(0,
                              nrow = ncol(X), ncol = p),
                   Precision = array(diag(0.01, ncol(X)),
                                     dim = c(ncol(X), ncol(X), p)))

library(abind)
result_VI2 <- get_CAVI(Y = Y, 
                       X = X,
                       q = params_Jags$q,
                       seed = 123, 
                       n_steps = 30, 
                       debug = FALSE, 
                       updates = c(Lambda = TRUE, Sigma = TRUE,
                                   Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                   Beta = TRUE, Z = TRUE),
                       get_ELBO_freq = 10,
                       params=params_Jags,
                       priors=priors
)
plot(result_VI2$ELBOS[-c(1:2),])
get_ELBO(Y=Y, params=params_Jags, priors=priors, X = X, XprimeX = t(X) %*% X)
get_ELBO(Y=Y, params=result_VI$params, priors=priors, X = X, XprimeX = t(X) %*% X)
get_ELBO(Y=Y, params=result_VI2$params, priors=priors, X = X, XprimeX = t(X) %*% X)

# Graphiques comparatifs ---------------------------------------------------------------

# test_tab <- tibble(groupe = rep(c("CAVI", "JAGS"), 2),
#                    params = rep(c("B1", "B2"), each = 2),
#                    mean = rnorm(4),
#                    sd = abs(rnorm(4)))
# test_tab
# ggplot() +
#   geom_function(fun = function(x) dnorm(x, mean = 2, sd = 4), colour = "red")

all_z_posterior <- expand.grid(i = 1:params_Jags$n, j = 1:params_Jags$p) %>% 
  pmap_dfr(function(i, j){
    tibble(method = c("VI", "Jags"),
           params = paste0("Z[", i, ",", j, "]"),
           mean = c(result_VI$params$Z$M[i, j],
                    params_Jags$Z$M[i, j]),
           sd = sqrt(c(result_VI$params$Z$S2[i, j],
                       params_Jags$Z$S2[i, j])))
  })

layers_list <- filter(all_z_posterior,
                      str_detect(params, "Z\\[1,")) %>%
  pmap(function(method, params, mean, sd) {
    geom_function(data = mutate(all_z_posterior, params = .env$params,
                                groupe = .env$method),
                  fun = dnorm,
                  mapping = aes(color = groupe),
                  args = list(mean = mean, sd = sd)
    )
  })
ggplot() +
  layers_list +
  facet_wrap(~params) +
  xlim(-5, 5)


all_beta_posterior <- expand.grid(i = 1:params_Jags$F_x, 
                                  j = 1:params_Jags$p) %>% 
  pmap_dfr(function(i, j){
    tibble(method = c("VI", "Jags"),
           params = paste0("Beta[", i, ",", j, "]"),
           mean = c(result_VI$params$Beta$M[i, j],
                    params_Jags$Beta$M[i, j]),
           sd = sqrt(c(result_VI$params$Beta$Cov[i, i, j],
                       params_Jags$Beta$Cov[i, i, j])))
  })

layers_list <- filter(all_beta_posterior) %>%
  pmap(function(method, params, mean, sd) {
    geom_function(data = mutate(all_beta_posterior, params = .env$params,
                                groupe = .env$method),
                  fun = dnorm,
                  mapping = aes(color = groupe),
                  args = list(mean = mean, sd = sd)
    )
  })
ggplot() +
  layers_list +
  facet_wrap(~params, nrow = 10) +
  xlim(-5, 5) 

all_sigmas_posterior <- map_dfr(1:params_Jags$p, function(j){
    tibble(method = c("VI", "Jags"),
           params = paste0("sigmas_m2[", j, "]"),
           shape = c(result_VI$params$Sigma$A[j],
                     params_Jags$Sigma$A[j]),
           rate = c(result_VI$params$Sigma$B[j],
                    params_Jags$Sigma$B[j]))
  })

layers_list <- filter(all_sigmas_posterior) %>%
  pmap(function(method, params, shape, rate) {
    geom_function(data = mutate(all_sigmas_posterior, params = .env$params,
                                groupe = .env$method),
                  fun = dgamma,
                  mapping = aes(color = groupe),
                  args = list(shape = shape, rate = rate)
    )
  })
ggplot() +
  layers_list +
  facet_wrap(~params, nrow = 10) +
  xlim(0, 10) 

1 / true_params$sigma2s
(params_Jags$Sigma$A / params_Jags$Sigma$B)  %>% round(2)
(result_VI$params$Sigma$A / result_VI$params$Sigma$B) %>% round(2)
1/true_params$sigma2s

(t(params_Jags$Lambda$M) %*% params_Jags$Lambda$M) %>% round(2)
(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M) %>% round(2)
(t(params_Jags$Eta$M) %*% params_Jags$Lambda$M) %>% round(2) %>% head()
(t(result_VI$params$Eta$M) %*% result_VI$params$Lambda$M) %>% round(2) %>% head()
