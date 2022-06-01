
# Cleaning environment ----------------------------------------------------

rm(list = ls())


# Librairies --------------------------------------------------------------

library(tidyverse)
library(mixtools)

# Problem data ------------------------------------------------------------

n <- 1000
p <- 5
q <- 3

# Generating latent variables ---------------------------------------------

sigma_m2s_true <- seq(10, 100, length.out = p)
sigma2s_true <- 1 / sigma_m2s_true 
# Residuals ---------------------------------------------------------------

set.seed(123)

residus <- sapply(sigma2s_true, function(s2) rnorm(n, 0, sqrt(s2)))


# Latent variables --------------------------------------------------------

set.seed(1234)
Eta_true <- matrix(rnorm(n * q), nrow = n, ncol = q) %>% 
  scale(center = FALSE) %>% 
  {. %*% svd(.)$v}

# Loadings ----------------------------------------------------------------

set.seed(12345)
Lambda_true <- matrix(sample(-1:1, size = p * q, replace = TRUE),
                 nrow = p, ncol = q)

# Observations

Y <- Eta_true %*% t(Lambda_true) + residus
rm(residus, q, n, p)
write.table(Y, file = "data_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
