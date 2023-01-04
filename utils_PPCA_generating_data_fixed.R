
# Cleaning environment ----------------------------------------------------

rm(list = ls())


# Librairies --------------------------------------------------------------

library(tidyverse)
library(mixtools)

# Problem data ------------------------------------------------------------

n <- 100
p <- 10
q_true <- 3
qk <-2

# Generating latent variables ---------------------------------------------

sigma_m2s_true <- seq(10, 100, length.out = p)
sigma2s_true <- 1 / sigma_m2s_true 
# Residuals ---------------------------------------------------------------

set.seed(123)

residus <- sapply(sigma2s_true, function(s2) rnorm(n, 0, sqrt(s2)))


# Latent variables --------------------------------------------------------

set.seed(1234)
Eta_true <- matrix(sample(-1:1, size = n * q_true, replace = TRUE),
                   nrow = n, ncol = q_true) 
# %>% 
#   scale(center = FALSE) %>% 
#   {. %*% svd(.)$v}

# Loadings ----------------------------------------------------------------

set.seed(12345)
Lambda_true <- matrix(sample(-1:1, size = p * q_true, replace = TRUE),
                 nrow = p, ncol = q_true)

# Fixed effects----------------------------------------------------------------
X <- matrix(sample(seq(fr=-10, to=10,by=0.1), size = n * qk, replace = TRUE),
                   nrow = n, ncol = qk)
beta_true <- matrix(runif(n = p*qk),
                      nrow = p, ncol = qk)*
  sqrt(matrix(sigma2s_true,
              nrow = p, ncol = qk))*10
beta_true=round(beta_true, dig=1)

###according to prior with c=100

# Observations

Y <- X %*% t(beta_true)+Eta_true %*% t(Lambda_true) + residus
rm(residus, q_true, n, p,qk)
write.table(Y, file = "data_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
write.table(X, file = "fixed_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
