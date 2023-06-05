
rm(list = ls())
library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(pacman) # For progress bar
library(torch)
source("utils_generating_data.R") # Creates experiment_params.rds


# Data --------------------------------------------------------------------

Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
true_params <- readRDS("experiment_params.rds")

# VI inference functions -----------------------------------------------------

source("utils_Poisson_PPCA_VI_functions.R")

result_VI <- get_CAVI(Y = Y, 
                      X = X,
                      q = 7,
                      seed = 5, 
                      n_steps = 2000, 
                      batch_prop = .4,
                      debug = FALSE, 
                      amortize = TRUE,
                      updates = c(Lambda = TRUE, Sigma = TRUE,
                                  Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                  Beta = TRUE, Z = TRUE),
                      get_ELBO_freq = 10)
result_VI$ELBOS %>% tail()
plot(result_VI$ELBOS[-c(1:2),])

result_VI$params$Sigma$B
save(result_VI, file = "result_VI.RData")

Z_predict <- X %*% result_VI$params$Beta$M + 
  t(result_VI$params$Eta$M) %*% result_VI$params$Lambda$M
Z_predict %>% round(1) %>% head()
(X %*% true_params$beta + true_params$Eta %*% t(true_params$Lambda)) %>% 
  round(1) %>% 
  head()

# Beta ? ------------------------------------------------------------------

true_params$beta
round(result_VI$params$Beta$M, 2)
result_VI$params$Beta$Cov
# Lambda? -----------------------------------------------------------------

(true_params$Lambda %*% t(true_params$Lambda)) %>% .[1:5,1:5] 
(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M)%>% .[1:5,1:5]

true_params$Lambda %*% t(true_params$Lambda) %>% 
  cov2cor()%>% .[1:5,1:5]
(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M) %>% 
  cov2cor()%>% .[1:5,1:5]

# Eta x Lambda ? ----------------------------------------------------------

head(true_params$Eta %*% t(true_params$Lambda)) %>% 
  round(2)
round(t(result_VI$params$Eta$M) %*% result_VI$params$Lambda$M, 2) %>% 
  head()

(round(true_params$Lambda %*% t(true_params$Lambda), 3) +
    diag(true_params$sigma2s)) %>% 
  cov2cor()%>% .[1:5,1:5]
(round(t(result_VI$params$Lambda$M) %*% result_VI$params$Lambda$M, 3) +
    diag(result_VI$params$Sigma$B / result_VI$params$Sigma$A)) %>% 
  cov2cor()%>% .[1:5,1:5]

range(result_VI$params$Z$M - X %*% result_VI$params$Beta$M)

plot(result_VI$ELBOS[-c(1:2),])
result_VI$params$Lambda$M
result_VI$params$Lambda$M %>% apply(1,var) %>% plot()

((result_VI$params$Sigma$A)/(result_VI$params$Sigma$B))^-0.5
true_params$sigma2s
((result_VI$params$Sigma$A)/(result_VI$params$Sigma$B))^-0.5%>% plot()


all_results <- lapply(1:50,
                      FUN = function(i)
                        #    get_result(Y = Y, X = X_data, n_steps = 50, seed = i)
                        get_CAVI(Y = Y, 
                                 X = X,
                                 q = 7,
                                 seed = 5*i, 
                                 n_steps = 30, 
                                 seed = i, 
                                 n_steps = 60, 
                                 debug = FALSE)
)

map_dfr(all_results, "ELBOS", .id = "Replicate") %>%
  filter(iteration > 4) %>%
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line() #+
 # theme(legend.position = "none")

best_seeds <- map_dfr(all_results, "ELBOS", .id = "Replicate") %>% 
  filter(iteration == max(iteration)) %>% 
  arrange(desc(ELBO)) %>% 
  head() %>% 
  mutate(seed = 5 * as.numeric(Replicate)) %>% 
  pull(seed)

all_results <- lapply(best_seeds,
                      FUN = function(i)
                        #    get_result(Y = Y, X = X_data, n_steps = 50, seed = i)
                        get_CAVI(Y = Y, 
                                 X = X,
                                 q = 7,
                                 seed = i, 
                                 n_steps = 100, 
                                 debug = FALSE)
)

all_params <- map(all_results, "params")
map(all_params, "Eta") %>% 
  map("Cov") %>% 
  map(range)
