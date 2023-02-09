
# Cleaning environment ----------------------------------------------------

rm(list = ls())


# Librairies --------------------------------------------------------------

library(tidyverse)
library(mixtools)

# Problem data ------------------------------------------------------------

n <- 100
p <- 10
q_true <- 3
F_x <-2

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
Lambda_true <- matrix(sample(-1:1, size = p * q_true, replace = TRUE) *0.1,
                      nrow = p, ncol = q_true)

# Fixed effects----------------------------------------------------------------
X <- matrix(sample(seq(from = -3, to = 3, by = 0.1), size = n * F_x, replace = TRUE),
            nrow = n, ncol = F_x) 
beta_true <- matrix(sample(-1:1, size = p * F_x, replace = TRUE) * 0.5,
                    nrow = p, ncol = F_x)

###according to prior with c=100

# Observations

Z <- X %*% t(beta_true) + Eta_true %*% t(Lambda_true) + residus
set.seed(123)
Y <- apply(Z, 2, function(z) rpois(length(z), exp(z)))
max(Y)
sum(Y == 0)
rm(residus, q_true, n, p,F_x)
write.table(Y, file = "data_Poisson_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
write.table(X, file = "fixed_Poisson_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
