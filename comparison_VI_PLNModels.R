rm(list = ls())
# Librairies --------------------------------------------------------------

library(tidyverse) # Pour la manipulation des données
library(PLNmodels) # Pour le modèle poisson log normal
library(future) # Pour calcul sur plusieurs processeurs
library(tidyverse)
library(abind) # To gather results together
library(parallel)
source("utils_generating_data.R") # Creates experiment_params.rds
source("utils_Poisson_PPCA_VI_functions.R")

# Chargement donnees simulees -----------------------------------------------


Y <- read.table("data_sets/synthetic/data_Y_Poisson_PPCA_with_covariates.txt", sep = ";") %>%
  as.matrix()
X <- read.table("data_sets/synthetic/data_covariates_Poisson_PPCA.txt", sep = ";") %>%
  as.matrix()
true_params <- readRDS("experiment_params.rds")



# Transformation pour la librairie ----------------------------------------

pln_data <- prepare_data(counts = Y, covariates = X[,-1,drop = FALSE],
                         offset = "none")

# PLN PCA -----------------------------------------------------------------

# Le modèle ajusté ici est expliqué brievement la
# https://pln-team.github.io/PLNmodels/articles/PLNPCA.html pour la PCA
# https://pln-team.github.io/PLNmodels/articles/PLNnetwork.html pour la partie réseau

# plan(multisession, workers = 3)
PCA_models <- PLNPCA(
  Abundance ~ V2,
  data  = pln_data,
  ranks = 1:10
)

plot(PCA_models)
PLN_res <- PLNmodels::getBestModel(PCA_models, crit = "BIC")
plot(PLN_res)


# VI PPCA -----------------------------------------------------------------

result_VI <- get_CAVI(Y = Y, 
                   X = X, 
                   q = 10,
                   seed = 1,
                   n_steps = 50, 
                   debug = FALSE)

plot(result_VI$ELBOS[-c(1:2),])
result_VI$params$Lambda$M

apply(result_VI$params$Lambda$M,1,var)

all_results <- lapply(1:4,
                      FUN = function(i)
                        #    get_result(Y = Y, X = X_data, n_steps = 50, seed = i)
                        get_CAVI(Y=Y_VIPCA, 
                                 X = X_VIPCA, 
                                 q = 5,
                                 seed = i, 
                                 n_steps = 100, 
                                 debug = FALSE)
)

map_dfr(all_results, "ELBOS", .id = "Replicate") %>%
  filter(iteration > 40) %>%
  ggplot(aes(x = iteration, y = ELBO, color = Replicate)) +
  geom_line() +
  theme(legend.position = "none")
# 
# get_result <- function(Y, X, q, seed, n_steps, init_params = NULL,
#                        priors = NULL,
#                        updates = c(Lambda = TRUE, Sigma = TRUE,
#                                    Eta = TRUE, Delta = TRUE, 
#                                    Phi = TRUE, Beta = !is.null(X), Z = TRUE), 
#                        debug = FALSE){
#   set.seed(seed)
#   if(is.null(X)){
#     F_x = 1
#   }
#   else{
#     F_x <- ncol(X)
#   }
#   p <- ncol(Y)
#   n <- nrow(Y)
#   # if(is.null(priors)){
#   #   priors <- list(Sigma = list(A = 3, B = 2), 
#   #                  Phi = list(A = 3/2, B = 3/2),
#   #                  Delta= list(A = c(5, rep(2, q - 1)), 
#   #                              B = 1),
#   #                  Beta = list(M = rep(0, ncol(X)),
#   #                              C = rep(0.01, ncol(X))))
#   # }
#   if(is.null(init_params)){
#     init_params <- list(Lambda = list(M = matrix(rnorm(p * q),
#                                                  nrow = q, ncol = p),
#                                       Cov = array(diag(1, q), 
#                                                   dim = c(q, q, p))),
#                         Beta = list(M = matrix(rnorm(p * F_x),
#                                                nrow = F_x, ncol = p),
#                                     Cov = array(diag(1, F_x), 
#                                                 dim = c(F_x, F_x, p))),
#                         Eta = list(M = matrix(rnorm(n * q), 
#                                               nrow = q, ncol = n),
#                                    Cov = array(diag(1, q), 
#                                                dim = c(q, q, n))),
#                         Sigma = list(A = rep(priors$Sigma$A + n / 2, p),
#                                      B = runif(p, 1, 5)),
#                         Delta = list(A = priors$Delta$A + 0.5 * p * (q + 1 - (1:q)),
#                                      B = runif(q, 1, 10)),
#                         Phi = list(A =  matrix(priors$Phi$A + 0.5, p, q),
#                                    B = matrix(runif(p * q, 1, 3), p, q)),
#                         Z = list(M = log(Y + 1),
#                                  S2 = matrix(.1, nrow = nrow(Y),
#                                              ncol = ncol(Y)))) 
#   }
#   result <- get_CAVI(data_ = Y, X = X, q = q, 
#                      n_steps = n_steps, params = init_params,
#                      updates = updates,
#                      priors = priors, debug = debug)
#   result
# }
# result = get_result(Y=Y_VIPCA, 
#                     X = X_VIPCA, 
#                     q = 5,
#                     seed = 1, 
#                     priors = list(Sigma = list(A = 3, B = 2), 
#                                   Phi = list(A = 3/2, B = 3/2),
#                                   Delta= list(A = c(5, rep(2, 4)), 
#                                               B = 1),
#                                   Beta = list(M = rep(0, 5),
#                                               C = rep(0.01, 5))),
#                     n_steps = 50, 
#                     updates = c(Lambda = TRUE, Sigma = TRUE,
#                                 Eta = TRUE, Delta = TRUE, 
#                                 Phi = TRUE, Beta = TRUE, Z = TRUE),
#                     debug = FALSE)
# 
