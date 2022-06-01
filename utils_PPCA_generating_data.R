
# Cleaning environment ----------------------------------------------------

rm(list = ls())


# Librairies --------------------------------------------------------------

library(tidyverse)
library(mixtools)

# Problem data ------------------------------------------------------------

n <- 1000
p <- 10
k <- 4

# Generating latent variables ---------------------------------------------

sigma2s_true <- seq(0.1, 0.01, length.out = p)

# Residuals ---------------------------------------------------------------

set.seed(123)

residus <- sapply(sigma2s_true, function(s2) rnorm(n, 0, sqrt(s2)))


# Latent variables --------------------------------------------------------

set.seed(1234)
Eta_true <- matrix(rnorm(n * k), nrow = n, ncol = k) %>% 
  scale(center = FALSE) %>% 
  {. %*% svd(.)$v}

# Loadings ----------------------------------------------------------------

set.seed(12345)
Lambda_true <- matrix(sample(-1:1, size = p * k, replace = TRUE),
                 nrow = p, ncol = k)

# Observations

Y <- Eta_true %*% t(Lambda_true) + residus
rm(residus, k, n, p)
write.table(Y, file = "data_PPCA.txt", sep = ";", col.names = FALSE,
            row.names = FALSE)
