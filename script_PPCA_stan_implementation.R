rm(list = ls())
library(tidyverse)
library(mixtools)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(rstan)
source("utils_PPCA_generating_data.R") # For true values


# Data --------------------------------------------------------------------

Y <- read.table("data_PPCA.txt", sep = ";") %>%
  as.matrix()
n <- nrow(Y); p <- ncol(Y)
q_max <- ncol(Y)
stan_data <- list(n = n, p = p, q_max = q_max, Y = Y)
pca <- "
	data {
		int<lower=0> n; // Number of samples
		int<lower=0> p; // The original dimension
		int<lower=0> q_max; // The latent dimension
		matrix[n, p] Y; // The data matrix
	}

	parameters {
		matrix[q_max, n] Z; // The latent matrix
		matrix[p, q_max] Lambda; // The weight matrix
		matrix<lower = 0>[p, q_max] phi; // Global precision matrix
		vector<lower = 0>[q_max] delta;
		vector<lower = 0>[p] sigma_diag;
		real<lower = 2> a1;
		real<lower = 3> a2;
	}
	transformed parameters{
	// Declarations
		vector<lower=0>[q_max] tau;
		vector<lower=0>[q_max] inv_tau;
		matrix<lower = 0>[p, q_max] inv_phi;
	//	matrix<lower = 0>[p, p] Sigma;
	// Assignations
		tau = exp(cumulative_sum(log(delta)));
		for(h in 1:q_max){
		  inv_tau[h] = 1 / tau[h];
		  for(j in 1:p){
		    inv_phi[j, h] = 1 / phi[j, h]; 
		  }
		}
	//	Sigma = diag_matrix(sigma_diag);
	}
	model {
	  delta[1] ~ gamma(a1, 1);
	  for(h in 2:q_max){
	    delta[h] ~ gamma(a2, 1);
	  }
		to_vector(phi) ~ gamma(3 * 0.5, 3 * 0.5);
		to_vector(Z) ~ normal(0,1);
		sigma_diag ~ inv_gamma(1, 1);
		for(h in 1:q_max){
		  Lambda[, h] ~ normal(0, sqrt(inv_phi[, h] * inv_tau[h]));
		}
		for(i in 1:n){
		  Y[i, ] ~ normal((Lambda * Z[,i])', sigma_diag);
		}
	} "




m <- stan_model(model_code = pca)
stan.fit <- rstan::sampling(m, data = stan_data, chains = 1, iter = 1e4,
                            pars = c("Lambda", "tau", "delta", "a1", "a2", "sigma_diag"))
stan.fit.vb <- vb(m, data = stan_data, 
                  algorithm = "meanfield",
                  pars = c("Lambda", "tau", "delta", "a1", "a2", "sigma_diag"))
all_res <- ggs(stan.fit.vb) 
all_res_mcmc <- ggs(stan.fit)
levels_Lambdas <- paste0("Lambda.", rep(1:p, q_max), ".", rep(1:q_max, each = p))

all_Lambdas <- all_res %>%   
  filter(str_detect(Parameter, "Lambda")) %>% 
  mutate(Parameter = factor(Parameter, levels = levels_Lambdas))
Lambda_est <- all_Lambdas %>% 
  group_by(Parameter) %>% 
  summarise(Estimate = mean(value)) %>% 
  as.data.frame() %>% 
  pull(Estimate) %>% 
  matrix(nrow = p)
apply(Lambda_est, 2, var) %>% plot()

LL_T <- (t(Lambda_est) %*% Lambda_est) %>% 
  matrix(nrow = q_max, ncol = q_max, 
         dimnames = list(paste0("S_", 1:q_max), paste0("S_", 1:q_max))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "row") %>% 
  pivot_longer(-row, names_to = "col", values_to = "value") %>% 
  mutate(row = factor(row, levels = rev(paste0("S_", 1:q_max))),
         col = factor(col, levels = paste0("S_", 1:q_max))) %>% 
  pull(value) %>% 
  matrix(nrow = p, ncol = q_max, byrow = TRUE)
  ggplot(aes(y = row, x = col, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c()
