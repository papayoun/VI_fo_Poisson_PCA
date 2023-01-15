

get_update_VI_Lambda <- function(Y, params){
  n <- nrow(Y)
  p <- ncol(Y)
  Eta <- params$Eta
  Sigma <- params$Sigma
  Phi <- params$Phi
  Delta <- params$Delta
  E_eta_prime_eta <- apply(Eta$Cov, c(1, 2), sum) + # Sum of variances
    Reduce(f = "+", 
           x = map(1:n, function(i){
             Eta$M[, i] %*% t(Eta$M[, i])
           }))
  # Calcul des precisions
  V_Lambda <- lapply(1:p, function(j){
    precision <- Sigma$A[j] / Sigma$B[j] * E_eta_prime_eta + 
      diag(Phi$A[j, ] / Phi$B[j, ] * cumprod(Delta$A) / cumprod(Delta$B))
    variance <- solve(precision)
    return(variance)
  }) %>% 
    abind(along = 3) # Mise au format array
  M_Lambda <- sapply(1:p, function(j){
    Sigma$A[j] / Sigma$B[j]  * V_Lambda[,, j] %*% Eta$M %*% Y[, j]
  })
  list(M = M_Lambda, Cov = V_Lambda)
}


get_update_VI_Eta <- function(Y, params){
  # Useful quantities
  # browser()
  n <- nrow(Y)
  p <- ncol(Y)
  Lambda <- params$Lambda
  q <- nrow(Lambda$M)
  Sigma <- params$Sigma
  Phi <- params$Phi
  Delta <- params$Delta
  
  #First, get the common covariance matrix
  common_cov <- lapply(1:p, function(j){
    Sigma$A[j] / Sigma$B[j] * (Lambda$M[,j] %*% t(Lambda$M[,j]) + Lambda$Cov[,,j])
  }) %>% 
    Reduce(f = "+") %>% 
    {. + diag(q)} %>% 
    solve()
  # The common matrix for all means
  mean_matrix <- common_cov %*% Lambda$M %*% diag(Sigma$A/Sigma$B)
  # The means
  M <- sapply(1:n, function(i){
    mean_matrix %*% Y[i, ] 
  })
  Cov <- array(common_cov, dim = c(q, q, n))
  list(M = M, Cov = Cov)
}

get_update_VI_Sigma_without_fixed_effects <- function(Y, params, priors){
  n <- nrow(Y)
  p <- ncol(Y)
  Lambda <- params$Lambda
  Eta <- params$Eta
  Phi <- params$Phi
  Delta <- params$Delta
  A <- rep(priors$Sigma$A + n / 2, p)
  E_eta_prime_eta <- apply(Eta$Cov, c(1, 2), sum) + # Sum of variances
    Reduce(f = "+", 
           x = lapply(1:n, function(i){
             Eta$M[, i] %*% t(Eta$M[, i])
           }))
  get_B_sigma_j <- function(j){
    term1 <- sum(0.5 * sum(Y[, j] * Y[, j]))
    term2 <- -sum(Y[, j] * (t(Eta$M) %*% Lambda$M[,j])) # Eta$M is coded in q x n
    term3 <- 0.5 * sum(E_eta_prime_eta * 
                         (Lambda$Cov[,, j] + Lambda$M[, j] %*% t(Lambda$M[, j]))) 
    term1 + term2 + term3
  }
  B <- priors$Sigma$B + sapply(1:p, get_B_sigma_j)
  list(A = A, B = B)
}

get_update_VI_Sigma <- function(Y, params, priors,
                                      X = 0, XprimeX = 0){
  # Useful quantities
  n <- nrow(Y)
  p <- ncol(Y)
  F_x <- ncol(X)
  Lambda <- params$Lambda
  Eta <- params$Eta
  Beta <- params$Beta
  # Posterior variationel de A
  update_without_X <- get_update_VI_Sigma_without_fixed_effects(Y, params, priors)
  if(all(X == 0)){ # Case when no covariates
    return(list(A = update_without_X$A, B = update_without_X$B))
  }
  else{
    A <- update_without_X$A + F_x * 0.5
    get_updade_term_B_sigma_j <- function(j){
      term1 <- 0.5 * sum((XprimeX + priors$Beta$C) * (Beta$Cov[,,j] + Beta$M[,j] %*% t(Beta$M[,j])))
      term2 <- -sum((Y[, j] -  t(Eta$M) %*% Lambda$M[,j]) * (X %*% Beta$M[,j])) # Eta$M is coded in q x n
      term1 + term2
    }
    B <- update_without_X$B + map_dbl(1:p, get_updade_term_B_sigma_j) 
  }
  return(list(A = A, B = B))
}

get_update_VI_Beta <- function(Y, params, priors, X, XprimeX){
  n <- nrow(Y)
  p <- ncol(Y)
  Sigma <- params$Sigma
  # Calcul des variances
  V_Beta <- lapply(1:p, function(j){
    precision <- Sigma$A[j] / Sigma$B[j] * (XprimeX + priors$Beta$C)
    variance <- solve(precision)
    return(variance)
  }) %>% 
    abind(along = 3) # Mise au format array
  M_Beta <- sapply(1:p, function(j){
    Sigma$A[j] / Sigma$B[j]  * V_Beta[,, j] %*% t(X) %*% (Y[, j] - t(params$Eta$M) %*% params$Lambda$M[, j])
  })
  list(M = M_Beta, Cov = V_Beta)
}

get_update_VI_Phi <- function(params, priors){
  Lambda <- params$Lambda
  q <- nrow(Lambda$M)
  Eta <- params$Eta
  Delta <- params$Delta
  
  # le meme pour tous
  A <- matrix(priors$Phi$A + 0.5, p, q)
  
  cumprod_Delta  = cumprod(Delta$A/Delta$B)
  B <- sapply(1:q, function(h){
    priors$Phi$B + 
      .5 * (Lambda$M[h, ]^2 + Lambda$Cov[h, h, ]) * cumprod_Delta[h]
  })
  list(A = A, B = B)
}

get_update_VI_Delta <- function(params, priors){
  p <- ncol(params$Lambda$M)
  q <- nrow(params$Lambda$M)
  new_A <- priors$Delta$A + 0.5 * p * (q + 1 - (1:q))
  E_phi_L2 <- sapply(1:q, function(h){
    E_phi_h <- params$Phi$A[, h] / params$Phi$B[, h]
    E_L2 <- params$Lambda$M[h, ]^2 + params$Lambda$Cov[h, h, ]
    sum(E_phi_h * E_L2)
  })
  new_B <- params$Delta$B
  new_B[1] <- priors$Delta$B + .5 * sum(E_phi_L2 * cumprod(new_A / new_B) / (new_A[1] / new_B[1]))
  for(k in 2:q){
    E_delta_all <- cumprod(new_A / new_B)
    E_delta_k <- (new_A[k] / new_B[k])
    new_B[k] <- priors$Delta$B  + .5 * sum(E_phi_L2[k:q] * 
                                             E_delta_all[k:q] / E_delta_k )
  }
  list(A = new_A, B = new_B)
}

get_entropy_normal <- function(Cov){
  0.5 * log(det(Cov))
}

get_entropy_gamma <- function(A, B){
  A - log(B) + lgamma(A) + (1 - A) * digamma(A)
}

get_expectation_gamma <- function(A, B){
  A / B
}

get_log_expectation_gamma <- function(A, B){
  digamma(A) -  log(B)
}

get_ELBO <- function(Y, params, priors){
  # Variational entropy term
  variational_entropy <- sum(apply(params$Lambda$Cov, 3, 
                                   get_entropy_normal)) + # Lambda
    sum(apply(params$Eta$Cov, 3, get_entropy_normal)) + # Eta
    sum(get_entropy_gamma(params$Sigma$A, params$Sigma$B)) + # Sigma 
    sum(get_entropy_gamma(params$Delta$A, params$Delta$B)) + # Delta 
    sum(get_entropy_gamma(params$Phi$A, params$Phi$B)) # Phi
  if(all(X != 0)){ # Case with covariates X with parameter beta
    variational_entropy <- variational_entropy+
      sum(apply(params$Beta$Cov, 3, get_entropy_normal))  # Lambda
                  }
  # Usefull expectations
  expectations_log_sigma <- get_log_expectation_gamma(params$Sigma$A, params$Sigma$B)
  expectations_sigma <- get_expectation_gamma(params$Sigma$A, params$Sigma$B)
  expectations_log_delta <- get_log_expectation_gamma(params$Delta$A, params$Delta$B)
  expectations_delta <- get_expectation_gamma(params$Delta$A, params$Delta$B)
  expectations_log_phi <- get_log_expectation_gamma(params$Phi$A, params$Phi$B)
  expectations_phi <- get_expectation_gamma(params$Phi$A, params$Phi$B)
  # Likelihood term
  E_eta_prime_eta <- apply(params$Eta$Cov, c(1, 2), sum) + # Sum of variances
    Reduce(f = "+", 
           x = lapply(1:n, function(i){
             params$Eta$M[, i] %*% t(params$Eta$M[, i])
           }))
  get_E_quadr_form <- function(j){
    term1 <- sum(0.5 * sum(Y[, j] * Y[, j]))
    term2 <- - sum(Y[, j] * (t(params$Eta$M) %*% params$Lambda$M[,j])) # Eta$M is coded in q x n
    term3 <- 0.5 * sum(E_eta_prime_eta * 
                         (params$Lambda$Cov[,, j] + params$Lambda$M[, j] %*% t(params$Lambda$M[, j]))) 
    term2bis <-0
    term3bis <-0
    term4 <- 0
    if(all(X != 0)){ # Case with covariates X with parameter beta
      term2bis <- -sum(Y[, j] * (X %*% params$Beta$M[,j]) ) # Eta$M is coded in q x n
      term3bis <- 0.5 * sum(XprimeX * 
                           (params$Beta$Cov[,, j] + params$Beta$M[, j] %*% t(params$Beta$M[, j]))) 
      term4 <- sum((X %*% params$Beta$M[,j])*t(params$Lambda$M[,j]%*%params$Eta$M))
      }
    
    term1 + term2 + term3 +term2bis +term3bis+term4
  }
  likelihood_expectation <- sum(.5 * n * expectations_log_sigma -
                                  expectations_sigma * map_dbl(1:p, get_E_quadr_form))
  # Priors terms
  prior_sigma_expectation <- sum((priors$Sigma$A - 1) * expectations_log_sigma - 
                                   priors$Sigma$B * expectations_sigma)
  prior_phi_expectation <- sum((priors$Phi$A - 1) * expectations_log_phi - 
                                 priors$Phi$B * expectations_phi)
  prior_delta_expectation <- sum((priors$Delta$A - 1) * expectations_log_delta - 
                                   priors$Delta$B * expectations_delta)
  prior_eta_expectation <- -0.5 * sum(map_dbl(1:n, function(i){
    sum(params$Eta$M[,i]^2 + diag(params$Eta$Cov[,, i]))
  }))
  prior_lambda_expectation <- 0.5 * p * sum(cumsum(expectations_log_delta)) +
    0.5 * sum(expectations_log_phi) - 
    0.5 * sum(map_dbl(1:p, function(j){
      sum(cumprod(expectations_delta) * diag(expectations_phi[j, ]) * 
            (diag(params$Lambda$Cov[,, j]) + params$Lambda$M[, j]^2))
    }))
  # A retravailler XXXXXXXXXXXXXXXX C precision?ou Variance
  prior_beta_expectation<-0
  if(all(X != 0)){ 
  prior_beta_expectation <- 0.5 * F_x * sum(expectations_log_sigma) -
    0.5 * sum(map_dbl(1:p, function(j){
      expectations_sigma[j] *sum(diag(priors$Beta$C) * 
            (params$Beta$Cov[,, j] + diag(params$Beta$M[, j]^2)))
    }))
  }
  ELBO <- variational_entropy +
    likelihood_expectation +
    prior_delta_expectation +
    prior_eta_expectation +
    prior_lambda_expectation +
    prior_sigma_expectation +
    prior_phi_expectation +
    prior_beta_expectation
  return(ELBO)
}


get_CAVI <- function(data_, q, n_steps, 
                     X = NULL, 
                     set_seed = FALSE,
                     params = NULL,
                     updates = c(Lambda = TRUE, Sigma = TRUE,
                                 Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                                 Beta = TRUE),
                     priors = list(Sigma = list(A = 1, B = 3), 
                                   Phi = list(A = 3/2, B = 3/2),
                                   Delta= list(A = c(2, rep(3, q - 1)), 
                                               B = 1)),
                     debug = FALSE){
  p <- ncol(data_); n <- nrow(data_); 
  if(is.null(X)){
    updates["Beta"] <- FALSE
    X <- matrix(0, nrow = nrow(data_), ncol = 1)
  }
  F_x <- ncol(X)
  XprimeX <- t(X) %*% X
  if(set_seed){
    set.seed(set_seed)
  }
  if(is.null(params)){
    params <- list(Lambda = list(M = matrix(rnorm(q * p), q, p),
                                 Cov = array(diag(1, q), dim = c(q, q, p))),
                   Beta = list(M = matrix(rnorm(F_x * p),
                                          nrow = F_x, ncol = p),
                               Cov = array(diag(1, F_x), 
                                           dim = c(F_x, F_x, p))),
                   Eta = list(M = matrix(rnorm(q * n), q, n),
                              Cov = array(diag(1, q), dim = c(q, q, n))),
                   Sigma = list(A = runif(p, 1, 3),
                                B = runif(p, 1, 3)),
                   Delta = list(A = runif(q, 2, 5) ,
                                B = rep(1, q)),
                   Phi = list(A = matrix(3/2, p, q),
                              B = matrix(3/2, p, q)))
  }
  
  # Propagation
  # progess_bar <- txtProgressBar(min = 0, max = n_steps)
  current_ELBO <- get_ELBO(data_, params, priors)
  ELBOS <- data.frame(iteration = 0, 
                      ELBO = current_ELBO)
  print(X)
  print(params$Beta$M)
  for(step_ in 1:n_steps){
    # Lambdas
    if(updates["Lambda"]){
      params$Lambda <- get_update_VI_Lambda(data_ - X %*% params$Beta$M, params)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Lambda")
        }
        current_ELBO <- new_ELBO
      }
    }
    # Etas
    if(updates["Eta"]){
      params$Eta <- get_update_VI_Eta(data_ - X %*% params$Beta$M, params)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Eta")
        }
        current_ELBO <- new_ELBO
      }
    }
    # Sigma
    if(updates["Sigma"]){
      params$Sigma <- get_update_VI_Sigma(data_, params, priors,
                                          X = X, XprimeX = XprimeX)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Sigma")
        }
        current_ELBO <- new_ELBO
      }
    }
    if(updates["Beta"]){
      params$Beta <- get_update_VI_Beta(Y = data_, params, priors, 
                                        X = X, XprimeX = XprimeX)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Beta")
        }
        current_ELBO <- new_ELBO
      }
    }
    # Phis
    if(updates["Phi"]){
      params$Phi <- get_update_VI_Phi(params, priors)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Phi")
        }
        current_ELBO <- new_ELBO
      }
    }
    if(updates["Delta"]){
      # deltas
      params$Delta <- get_update_VI_Delta(params, priors)
      if(debug){
        new_ELBO <- get_ELBO(data_, params, priors)
        if(new_ELBO < current_ELBO){
          print(paste("Problem at iteration", step_, "after updating Delta"))
        }
        current_ELBO <- new_ELBO
      }
    }
     if(n_steps){
       ELBOS <- bind_rows(ELBOS,
                          data.frame(iteration = step_,
                                     ELBO = get_ELBO(data_, params, priors)))
     }
  }
  return(list(ELBOS = ELBOS, params = params))
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