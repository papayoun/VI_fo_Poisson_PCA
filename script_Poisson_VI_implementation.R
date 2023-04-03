rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(pacman) # For progress bar
# library(pacman) # For progress bar (ne marche pas... ou travail non fini par celui qui l'a ajouté!
# Du coup j'ai bricolé mon propre progress bar Eric)
source("utils_generating_data.R") # Creates experiment_params.rds


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
true_params <- readRDS("experiment_params.rds")

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_functions.R")

# p <- ncol(Y)
# q <- max(1, ncol(Y) - 1)
# n <- nrow(Y)
# priors = list(Sigma = list(A = 3, B = 2), 
#               Phi = list(A = 3/2, B = 3/2),
#               Delta= list(A = rep(3, q), 
#                           B = 1),
#               Beta = list(M = matrix(0,
#                                      nrow = ncol(X), ncol = p),
#                           Precision = array(diag(0.01, ncol(X)),
#                                             dim = c(ncol(X), ncol(X), p))))
result_VI <- get_CAVI(Y = Y, 
                      X = X,
                      q = 7,
                      seed = 123, 
                      n_steps = 30, 
                      debug = FALSE, 
                      updates = c(Lambda = TRUE, Sigma = TRUE,
                                  Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                  Beta = TRUE, Z = TRUE),
                      get_ELBO_freq = 10)
plot(result_VI$ELBOS[-c(1:2),])
save(result_VI, file = "result_VI.RData")

# Beta ? ------------------------------------------------------------------

true_params$beta
round(result_VI$params$Beta$M, 12)
result_VI$params$Beta$Cov
# Lambda? -----------------------------------------------------------------


true_params$Lambda %*% t(true_params$Lambda) %>% 
  cov2cor()
(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M) %>% 
  cov2cor()

# Eta x Lambda ? ----------------------------------------------------------

head(true_params$Eta %*% t(true_params$Lambda)) %>% 
  round(2)
round(t(result_VI$params$Eta$M) %*% result_VI$params$Lambda$M, 2) %>% 
  head()

(round(true_params$Lambda %*% t(true_params$Lambda), 3) +
    diag(true_params$sigma2s)) %>% 
  cov2cor()
(round(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M, 3) +
    diag(result_VI$params$Sigma$B / result_VI$params$Sigma$A)) %>% 
  cov2cor()

range(result_VI$params$Z$M - X %*% result_VI$params$Beta$M)

plot(result_VI$ELBOS[-c(1:2),])
result_VI$params$Lambda$M

all_results <- lapply(1:4,
                      FUN = function(i)
                        #    get_result(Y = Y, X = X_data, n_steps = 50, seed = i)
                        get_CAVI(Y = Y, 
                                 X = X,
                                 q = 7,
                                 seed = i, 
                                 n_steps = 150, 
                                 debug = FALSE)
)

map_dfr(all_results, "ELBOS", .id = "Replicate") %>%
  filter(iteration > 40) %>%
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line() +
  theme(legend.position = "none")

all_params <- map(all_results, "params")
map(all_params, "Eta") %>% 
  map("Cov") %>% 
  map(range)