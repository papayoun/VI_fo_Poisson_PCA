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
data_for_JAGS <- list(
  n = dim(Y)[1],
  p = dim(Y)[2],
  q = dim(Y)[2] - 1,
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
  variable.names = c("beta","lambda","preciE", "eta"),
  n.iter = n.iter,
  thin = thin
)

result_jags <- ggs(jsamples)
save(result_jags, file = "result_jags.RData")
