rm(list = ls())
library(tidyverse)
library(mixtools)
library(abind) # To gather results together
library(parallel) # For some parallel computations
source("utils_PPCA_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_PPCA.txt", sep = ";") %>%
  as.matrix()

# Gibbs inference functions -----------------------------------------------------

source("utils_PPCA_gibbs_function.R")

# Gibbs parameters

n_steps <- 100

# Oracle case, Eta_known, k_tilde = 2 ------------------------------------------------

res_known_eta <- get_gibbs_sample(data_ = Y, 
                                  n_steps = n_steps, k_tilde = 2,
                           Eta_true = Eta_true)

## Estimations of Sigmas

Sigma_hat_df <- format_matrix(res_known_eta$Sigma, "sigma", "^2") 
ggplot(Sigma_hat_df) + 
  aes(x = iteration, y = Estimate) +
  geom_line() +
  facet_wrap(~Parameter, nrow = 2, labeller = label_parsed) +
  geom_hline(data = tibble(Parameter = unique(Sigma_hat_df$Parameter),
                           Truth = sigma2s_true),
             aes(yintercept = Truth),
             color = "red") +
  scale_y_continuous(trans = "log10")

## Estimations of Lambdas

Lambda_hat_df <- format_array(res_known_eta$Lambda, "Lambda") 
ggplot(Lambda_hat_df) + 
  aes(x = Estimate) +
  geom_density() +
  facet_wrap(~Parameter, nrow = nrow(Lambda_true), labeller = label_parsed) +
  geom_vline(data = tibble(Parameter = unique(Lambda_hat_df$Parameter),
                           Truth = as.numeric(t(Lambda_true))),
             aes(xintercept = Truth),
             color = "red")


# Oracle case, Lambda known, k = 2 ----------------------------------------

res_known_lambda <- get_gibbs_sample(data_ = Y, n_steps = n_steps,
                                     k_tilde = 3,
                           Lambda_true = Lambda_true)

## Estimation of Sigmas

Sigma_hat_df <- format_matrix(res_known_lambda$Sigma, "sigma", "^2") 
ggplot(Sigma_hat_df) + 
  aes(x = iteration, y = Estimate) +
  geom_line() +
  facet_wrap(~Parameter, nrow = 2, labeller = label_parsed) +
  geom_hline(data = tibble(Parameter = unique(Sigma_hat_df$Parameter),
                           Truth = sigma2s_true),
             aes(yintercept = Truth),
             color = "red") +
  scale_y_continuous(trans = "log10")

## Estimaton of Etas

### Choosing 10 rows of etas randomly

select_eta_rows <- sample(1:nrow(Y), size = 10, 
                          replace = FALSE) %>% 
  sort()
select_eta_rows <- 1:6

Eta_hat_df <- format_array(res_known_lambda$Eta, "eta", 
                           row_indexes = select_eta_rows) 
ggplot(Eta_hat_df) + 
  aes(x = Estimate) +
  geom_density() +
  facet_wrap(~Parameter, nrow = length(select_eta_rows),
             labeller = label_parsed) +
  geom_vline(data = tibble(Parameter = unique(Eta_hat_df$Parameter),
                           Truth = as.numeric(t(Eta_true[select_eta_rows, ]))),
             aes(xintercept = Truth),
             color = "red")


# All unknown case, k = 5 -------------------------------------------------

n_steps <- 10000
n_burn <- 50
n_thin <- 50

k_tilde <- 4

rm(res_known_eta, res_known_lambda, Sigma_hat_df, Lambda_hat_df, Eta_hat_df,
   select_eta_rows) 
result <- get_gibbs_sample(data_ = Y, n_steps = n_steps, k_tilde = k_tilde, 
                           burn = n_burn, thin = n_thin)

## Estmation of Lambda Lambda'


n_samples <- ncol(result$Sigma)
LL_t <- mclapply(1:n_samples, function(i){
  result$Lambda[,,i] %*% t(result$Lambda[,,i]) 
}) %>% 
  abind(along = 3)

LL_t_hat_df <- format_array(LL_t, "LL_t") 
LL_t_true <- Lambda_true %*% t(Lambda_true)

LL_t_hat_df %>% 
  ggplot(aes(x = iteration, y = Estimate, color = Parameter)) +
  geom_line()

ggplot(LL_t_hat_df) + 
  aes(x = Estimate) +
  geom_density() +
  facet_wrap(~Parameter, nrow =  nrow(Lambda_true), labeller = label_parsed,
             scales = "free_y") +
  geom_vline(data = tibble(Parameter = unique(LL_t_hat_df$Parameter),
                           Truth = as.numeric(t(LL_t_true))),
             aes(xintercept = Truth),
             color = "red")

## Estmation of Lambda

Lambda_hat_df <- format_array(result$Lambda,
                              "Lambda") 
Lambda_hat_df %>% 
  separate(Parameter, into = c("species", "feature"), sep = "-") %>% 
  mutate(feature = paste0("feature_", str_remove(feature , "]")),
         species = str_replace(species, "Lambda\\[", "species_")) %>% 
  ggplot(aes(x = iteration, y = Estimate, color = species)) +
  facet_wrap(~feature, scales = "free_y") +
  geom_line()

if(k_tilde > ncol(Lambda_true)){
  Lambda_true <- cbind(Lambda_true,
                       matrix(0, 
                              ncol = k_tilde - ncol(Lambda_true),
                              nrow = nrow(Lambda_true)))
}

ggplot(Lambda_hat_df) + 
  aes(x = Estimate) +
  geom_density() +
  facet_wrap(~Parameter, nrow =  nrow(Lambda_true), labeller = label_parsed,
             scales = "free_y") +
  geom_vline(data = tibble(Parameter = unique(Lambda_hat_df$Parameter),
                           Truth = as.numeric(t(Lambda_true))),
             aes(xintercept = Truth),
             color = "red")

## Representation of the loadings (species of niche)

Lambda_hat_df %>% 
  filter((str_detect(Parameter, "-1]") | str_detect(Parameter, "-2]"))) %>% 
  separate(Parameter, into = c("species", "feature"), sep = "-") %>% 
  mutate(feature = paste0("feature_", str_remove(feature , "]")),
         species = str_replace(species, "Lambda\\[", "species_")) %>%
  pivot_wider(names_from = "feature", values_from = "Estimate") %>% 
  ggplot(aes(x = feature_1, y = feature_2, color = species)) +
  geom_point() 

