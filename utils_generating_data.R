
# Cleaning environment ----------------------------------------------------

rm(list = ls())


# Librairies --------------------------------------------------------------

library(tidyverse)

# Problem data ------------------------------------------------------------

n <- 2000
p <- 10
q_true <- 2
F_x <- 5

# Generating latent variables ---------------------------------------------

sigma_m2s_true <- seq(5, 10, length.out = p)
sigma2s_true <- 1 / sigma_m2s_true 

# Residuals ---------------------------------------------------------------

set.seed(123)

residus <- sapply(sigma2s_true, function(s2) rnorm(n, 0, sqrt(s2)))

# Latent variables --------------------------------------------------------

set.seed(1234)
Eta_true <- matrix(sample(-1:1, size = n * q_true, replace = TRUE),
                   nrow = n, ncol = q_true) 
Eta_true <- matrix(rnorm(n * q_true),
                   nrow = n, ncol = q_true) 
# Loadings ----------------------------------------------------------------

set.seed(123)
Lambda_true <- matrix(sample(-1:1, size = p * q_true, replace = TRUE),
                 nrow = p, ncol = q_true)
Lambda_true <- matrix(round(rnorm(p * q_true), 3),
                   nrow = p, ncol = q_true) %>% 
  svd() %>% 
  {.$u * 2}
# Fixed effects----------------------------------------------------------------

X <- cbind(1,
           matrix(round(rnorm(n * (F_x - 1)), 2),
                  nrow = n, ncol = F_x - 1))

beta_true <- 0.3 * matrix(sample(-1:1, size = p * F_x, replace = TRUE),
                          nrow = F_x, ncol = p) +
  rbind(rep(2, p),
        matrix(0, nrow = F_x - 1, ncol = p))
# Normal PPCA -------------------------------------------------------------

# Without fixed

Y <- Eta_true %*% t(Lambda_true) + residus
write.table(round(Y, 4), file = "data_sets/synthetic/data_Y_Normal_PPCA.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)

# With fixed

write.table(round(Y + X %*% beta_true, 4), file = "data_sets/synthetic/data_Y_Normal_PPCA_with_covariates.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)

write.table(X, file = "data_sets/synthetic/data_covariates_Normal_PPCA.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)


# Poisson PPCA ------------------------------------------------------------

# Without fixed

Y <- apply(Eta_true %*% t(Lambda_true) + residus, 2,
           function(z) rpois(length(z), exp(z)))
write.table(Y, file = "data_sets/synthetic/data_Y_Poisson_PPCA.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)

# With fixed

Y <- apply(X %*% beta_true + Eta_true %*% t(Lambda_true) + residus, 2,
           function(z) rpois(length(z), exp(z)))
write.table(Y, file = "data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)
write.table(X, file = "data_sets/synthetic/data_covariates_Poisson_PPCA.txt", 
            sep = ";", col.names = FALSE,
            row.names = FALSE)

true_params <- list(beta = beta_true,
                    Eta = Eta_true,
                    Lambda = Lambda_true,
                    sigma2s = sigma2s_true)
saveRDS(true_params, file = "experiment_params.rds")
rm(list = ls())

