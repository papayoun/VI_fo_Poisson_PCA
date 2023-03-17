rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
# library(pacman) # For progress bar (ne marche pas... ou travail non fini par celui qui l'a ajouté!
# Du coup j'ai bricolé mon propre progress bar Eric)

# source("utils_Poisson_PPCA_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_functions.R")

result <- get_CAVI(Y = Y, 
                   X = X,
                   q = 7,
                   seed = 1, 
                   n_steps = 50, 
                   debug = FALSE)


plot(result$ELBOS[-c(1:2),])
result$params$Lambda$M
result_VI <- result
save(result_VI, file = "result_VI.RData")
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

# params_list <- map(all_results, "params")
# 
# best <- params_list[[1]]
# t(best$Beta$M) %>% 
#   round(2)
# beta_true
# t(best$Lambda$M)
# Lambda_true

