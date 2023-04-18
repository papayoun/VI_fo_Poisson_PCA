rm(list = ls())
library(ggmcmc)
library(rjags)
source("utils_generating_data.R") # Creates experiment_params.rds

# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", 
                sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", 
                sep = ";") %>%
  as.matrix()

#Un peu long mais faisable en jags pour n=100,p=10,q=7,F_x=2.
modelString = "
model{
#a1 <- 2
#a2 <- 3
a[1] ~ dgamma(2,1)T(2,)
a[2] ~ dgamma(2,1)T(3,)
deminu <- 3/2
a_sigma <- 0.1
b_sigma <- 0.1
precibeta <- 0.01

#shrinkage prior
for( h in 2:q){delta[h] ~ dgamma(a[2],1)}
              delta[1] ~ dgamma(a[1],1)
for( h in 1:q){ tau[h] <- prod(delta[1:h])} 
for( j in 1:p){
               for( h in 1:q){
                              phi[j,h] ~ dgamma(deminu,deminu)
                              precilambda[j,h] <- phi[j,h]*tau[h]
                              lambda[j,h]~ dnorm(0,precilambda[j,h])
                              }
               }
# N(0,1) prior for eta
for(i in 1:n){
              for (h in 1:q){
                      eta[i,h] ~ dnorm(0,1)
                            }
}
for( f in 1:F_x){
               for( j in 1:p){
               beta[j,f]~ dnorm(0,precibeta)
               }
}
muZ <- eta %*% t(lambda)+X %*% t(beta)

# latent Z
for(j in 1:p){ preciE[j] ~ dgamma(a_sigma,b_sigma)
              for (i in 1:n){
                      Z[i,j] ~ dnorm(muZ[i,j],preciE[j])
                            }
}
#Finally Y
for(j in 1:p){ 
              for (i in 1:n){ muY[i,j]<- exp(Z[i,j])
                      Y[i,j] ~ dpois(muY[i,j])
                            }

}
}"
#  q = dim(Y)[2] - 1
q <- 7 #(params$q)
data_for_JAGS <- list(
  n = dim(Y)[1],
  p = dim(Y)[2],
  q = q,
  Y = Y,
  X = X,
  F_x = ncol(X)
)

n.adapt = 300
burnin = 300
n.iter = burnin * 3
thin = 1
jm <- jags.model(
  file = textConnection(modelString),
  data = data_for_JAGS,
  n.chains = 3 ,
  n.adapt = n.adapt
)
update(jm, burnin)
jsamples <- coda.samples(
  model = jm,
  variable.names = c("beta","lambda","preciE", "eta","Z","delta","phi"),
  n.iter = n.iter,
  thin = thin
)

result_jags <- ggs(jsamples)
save(result_jags, file = "result_jags.RData")
load("result_jags.RData")
true_params <- readRDS("experiment_params.rds")
load("result_VI.RData")
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
params_Jags =list()
# Z
nbi=dim(result_VI$params$Z$M)[1]
nbj=dim(result_VI$params$Z$M)[2]
library(tidyverse)
result_jags %>% dplyr::filter(str_detect(Parameter,"Z")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  dplyr::select(c(-Iteration,-Chain)) ->tabZ
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
params_Jags$Phi$A<- MM^2/VV
params_Jags$Phi$B<- MM/VV
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

source("utils_Poisson_PPCA_VI_functions.R")

library(abind)
result_VI2 <- get_CAVI(Y = Y, 
                      X = X,
                      q = 7,
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
