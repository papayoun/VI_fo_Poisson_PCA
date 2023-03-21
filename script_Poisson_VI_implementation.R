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

result <- get_CAVI(Y = Y, 
                   X = X,
                   q = max(ncol(Y) - 1, 1),
                   seed = 1, 
                   n_steps = 50, 
                   debug = FALSE, 
                   updates = c(Lambda = TRUE, Sigma = TRUE,
                               Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                               Beta = TRUE, Z = TRUE),
                   get_ELBO_freq = 10)
result_VI <- result
save(result_VI, file = "result_VI.RData")

# Beta ? ------------------------------------------------------------------

true_params$beta
round(result$params$Beta$M, 12)

# Lambda? -----------------------------------------------------------------


true_params$Lambda %*% t(true_params$Lambda)
round(t(result$params$Lambda$M) %*% result$params$Lambda$M, 6)


# Eta x Lambda ? ----------------------------------------------------------

head(true_params$Eta %*% t(true_params$Lambda))
round(t(result$params$Eta$M) %*% result$params$Lambda$M, 2) %>% 
  head()

range(result$params$Z$M - X %*% result$params$Beta$M)

plot(result$ELBOS[-c(1:2),])
result$params$Lambda$M

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
