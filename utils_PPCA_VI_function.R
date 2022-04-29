get_update_VI_Lambda <- function(Y, Phi_mat, tau_vec, Eta_mat, sigma2_vec){
  sapply(1:length(sigma2_vec), function(j){
    Dj_inv <- diag(Phi_mat[j,] * tau_vec)
    cov_matrix <- solve(Dj_inv +
                          t(Eta_mat) %*% Eta_mat / sigma2_vec[j])
    mu <- cov_matrix %*% t(Eta_mat) %*% Y[, j] /  sigma2_vec[j]
    rmvnorm(n = 1, mu = as.numeric(mu), sigma = cov_matrix) %>%
      as.numeric()
  }) %>%
    t()
}

get_update_Sigma <- function(Y, Eta_mat, Lambda_mat, a_sigma, b_sigma){
  p <- nrow(Lambda_mat)
  1 / rgamma(p,
             a_sigma + .5 * nrow(Y),
             b_sigma + .5 * colSums((Y - Eta_mat %*% t(Lambda_mat))^2))
}

get_update_Eta <- function(Y, Lambda_mat, sigma2_vec){
  k_tilde <- ncol(Lambda_mat)
  cov_matrix_Etas <- solve(diag(1, k_tilde) +
                             t(Lambda_mat) %*% diag(1 / sigma2_vec) %*% Lambda_mat)
  mu_Etas <- cov_matrix_Etas %*% t(Lambda_mat) %*% diag(1 / sigma2_vec) %*% t(Y)
  sapply(1:nrow(Y), function(i){
    rmvnorm(n = 1, mu = mu_Etas[, i], sigma = cov_matrix_Etas) %>%
      as.numeric()
  }) %>%
    t()
}

get_update_Phi <- function(Lambda_mat, tau_vec, nu){
  p <- nrow(Lambda_mat)
  k_tilde <- ncol(Lambda_mat)
  rgamma(p * k_tilde, 
         .5 * (nu + 1), 
         .5 * (nu + Lambda_mat^2 * tau_vec))
}

get_tau_minus <- function(delta_vec, h){
  if(h == 1){
    cumprod(c(1, delta_vec[-1]))
  }
  else{
    (cumprod(delta_vec) / delta_vec[h])[-c(1:(h - 1))]
  }
}
get_update_deltas <- function(delta_vec, Lambda_mat, Phi_mat, a_1, a_2){
  output <- delta_vec
  Lambda_Phi <- colSums(Lambda_mat^2 * Phi_mat)
  p <- nrow(Lambda_mat)
  k_tilde <- ncol(Lambda_mat)
  output[1] <- rgamma(1, 
                      a_1 + .5 * p * k_tilde,
                      1 + .5 * sum((get_tau_minus(output, 1) * Lambda_Phi))) 
  for(h in 2:k_tilde){
    output[h] <- rgamma(1, 
                        a_2 + .5 * p * k_tilde,
                        1 + .5 * sum(get_tau_minus(output, h) * Lambda_Phi[-(1:(h - 1))]))
  }
  output
}

get_gibbs_sample <- function(data_, n_steps, 
                             k_tilde, 
                             raw_output = TRUE,
                             burn = 0,
                             thin = 1,
                             Lambda_true = NULL,
                             Eta_true = NULL,
                             sigma2s_true = NULL){
  
  # Needed quantities related to data
  
  p <- ncol(data_); n <- nrow(data_)
  
  # Hyperparameters ---------------------------------------------------------
  
  nu <- 3; a_sigma <- 3; b_sigma <- 2
 # a_1 <- 5; a_2 <- 2
  a_1 <- 3; a_2 <- 2
  # Creating outputs 
  
  # First, Lambdas
  if(!is.null(Lambda_true)){
    Lambdas <- array(Lambda_true, dim = c(p, k_tilde, n_steps + 1))
    update_Lambda = FALSE
  }
  else{
    Lambdas <- array(dim = c(p, k_tilde, n_steps + 1))
    update_Lambda = TRUE
  } 
  # Then, sigma2s
  if(!is.null(sigma2s_true)){
    sigma2s <- matrix(sigma2s_true, nrow = p, ncol = n_steps + 1)
    update_Sigma = FALSE
  }
  else{
    sigma2s <- matrix(nrow = p, ncol = n_steps + 1)
    update_Sigma = TRUE
  } 
  # Then, Etas
  if(!is.null(Eta_true)){
    Etas <- array(Eta_true, dim = c(n, k_tilde, n_steps + 1))
    update_Eta = FALSE
  }
  else{
    Etas <- array(dim = c(n, k_tilde, n_steps + 1))
    update_Eta = TRUE
  }  
  Phis <- array(dim = c(p, k_tilde, n_steps + 1))
  deltas <- matrix(nrow = k_tilde, ncol = n_steps + 1)
  taus <- matrix(nrow = k_tilde, ncol = n_steps + 1)
  
  # Initialisation of the chain
  # deltas
  deltas[1, 1] <- rgamma(1, a_1, 1)
  deltas[-1, 1] <- rgamma(k_tilde - 1, a_2, 1)
  # taus
  taus[, 1] <- cumprod(deltas[, 1])
  # Phis
  Phis[,, 1] <- rgamma(p * k_tilde, nu / 2, nu / 2)
  # sigma2s
  if(update_Sigma){
    sigma2s[, 1] <- 1 / sort(rgamma(p, a_sigma, b_sigma))
  }
  # Etas
  if(update_Eta){
    Etas[,, 1] <- rnorm(n * k_tilde)
  }
  # Lambda
  if(update_Lambda){
    Lambdas[,, 1] <- rnorm(p * k_tilde)
  }
  # Propagation
  progess_bar <- txtProgressBar(min = 0, max = n_steps)
  for(step in 1:n_steps){
    setTxtProgressBar(progess_bar, step)
    # Lambdas
    if(update_Lambda){
      Lambdas[,, step + 1] <- get_update_Lambda(Y = data_, 
                                                Phi_mat = Phis[,, step], 
                                                tau_vec = taus[, step], 
                                                Eta_mat = Etas[,, step], 
                                                sigma2_vec = sigma2s[, step])
    }
    # Sigma
    if(update_Sigma){
      sigma2s[, step + 1] <- get_update_Sigma(Y = data_,
                                              Eta_mat = Etas[,, step], 
                                              Lambda_mat = Lambdas[,, step + 1],
                                              a_sigma = a_sigma, b_sigma = b_sigma)
    }
    # Etas
    if(update_Eta){
      Etas[,, step + 1] <- get_update_Eta(Y = data_, 
                                          Lambda_mat = Lambdas[,, step + 1],
                                          sigma2_vec = sigma2s[, step + 1])
    }
    # Phis
    Phis[,, step + 1] <- get_update_Phi(Lambda_mat = Lambdas[,, step + 1], 
                                        tau_vec = taus[, step], 
                                        nu = nu)
    # deltas
    deltas[, step + 1] <- get_update_deltas(delta_vec = deltas[, step], 
                                            Lambda_mat = Lambdas[,, step + 1], 
                                            Phi_mat =  Phis[,, step + 1], 
                                            a_1 = a_1, a_2 = a_2)
    taus[, step + 1] <- cumprod(deltas[, step + 1])
  }
  # Thin and burn
  selection <- seq(from = burn + 1, to = n_steps + 1, by = thin)
  Lambdas <- Lambdas[,, selection]
  Etas <- Etas[,, selection]
  sigma2s <- sigma2s[, selection]
  deltas <- deltas[, selection]
  taus <- taus[, selection]
  Phis <- Phis[,, selection]
  # Return
  if(raw_output) 
    result = list(Lambda = Lambdas, Eta = Etas, Sigma = sigma2s,
                  deltas = deltas,
                  taus = taus, Phi = Phis)
  else{
    # Formatting results
    # Arrays
    Lambda_df <- format_array(Lambdas, "Lambda")
    Eta_df <- format_array(Etas, "eta")
    Phi_df <- format_array(Phis, "phi")
    # Matrices
    deltas_df <- format_matrix(deltas, "delta")
    taus_df <- format_matrix(taus, "tau")
    sigma2_df <- format_matrix(sigma2s, "sigma", "^2")
    result = bind_row(Lambda_df, Eta_df, Phi_df, deltas_df, taus_df, sigma2_df)
  }
  close(progess_bar)
  return(result)
}

# Function to format an array of estimated parameters to a data.frame
# easily handle by ggplot2

format_array <- function(array_, param_name, row_indexes = NULL){
  if(is.null(row_indexes)){
    row_indexes <- 1:(dim(array_)[1])
  }
  else{
    row_indexes <- sort(row_indexes)
    array_ <- array_[row_indexes,,]
  }
  dim1 <- length(row_indexes)
  dim2 <- dim(array_)[2]
  n_iteration <- dim(array_)[3]
  names_levels <- paste0(param_name, "[",
                         rep(row_indexes, each = dim2),
                         "-",
                         rep(1:dim2, dim1),
                         "]")
  map_dfr(1:n_iteration,
          function(l){
            array_[,, l] %>% 
              as.data.frame() %>% 
              setNames(paste0("col", 1:dim2)) %>% 
              mutate(iteration = l,
                     row = row_indexes) %>% 
              pivot_longer(-c("iteration", "row"), 
                           names_to = "column",
                           names_prefix = "col",
                           values_to = "Estimate") %>% 
              unite(Parameter, row, column, sep = '-') %>% 
              mutate(Parameter = paste0(param_name, "[", Parameter, "]"),
                     Parameter = factor(Parameter,
                                        levels = names_levels)) %>% 
              arrange(Parameter)
          })
}

# Function to format a matrix of estimated parameters to a data.frame
# easily handle by ggplot2
format_matrix <- function(matrix_, # Matrix of values d x n_iterations
                          param_name, # String giving parameter name
                          suffix_ = NULL){ # suffix of the name typically "^2" when
  # the name is sigma
  dim1 <- dim(matrix_)[1] # DImension of the parameter
  n_iteration <- dim(matrix_)[2]
  names_levels <- paste0(param_name, "[", 1:dim1, "]", suffix_)
  matrix_ %>% 
    as.data.frame() %>% 
    setNames(paste0("it", 1:n_iteration)) %>%
    mutate(Parameter = names_levels) %>%  
    pivot_longer(-Parameter,
                 names_to = "iteration",
                 names_prefix = "it", 
                 names_transform = list(iteration = as.numeric),
                 values_to = "Estimate") %>% 
    mutate(Parameter = factor(Parameter, level = names_levels)) %>% 
    arrange(Parameter)
}