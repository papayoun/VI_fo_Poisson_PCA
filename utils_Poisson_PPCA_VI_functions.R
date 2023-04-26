
# Update functions --------------------------------------------------------


get_update_Poisson_VI_Lambda <- function(params, X){
  Eta <- params$Eta
  Sigma <- params$Sigma
  Phi <- params$Phi
  Delta <- params$Delta
  E_eta_prime_eta <- apply(Eta$Cov, c(1, 2), sum) + # Sum of variances
    Reduce(f = "+", 
           x = map(1:params$n, function(i){
             Eta$M[, i] %*% t(Eta$M[, i])
           }))
  # Calcul des precisions
  V_Lambda <- lapply(1:params$p, function(j){
    precision <- Sigma$A[j] / Sigma$B[j] * E_eta_prime_eta + 
      diag(x = Phi$A[j, ] / Phi$B[j, ] * cumprod(Delta$A) / cumprod(Delta$B),
           nrow = params$q, ncol = params$q)
    variance <- solve(precision)
    return(variance)
  }) %>% 
    abind(along = 3) # Mise au format array
  M_Lambda <- sapply(1:params$p, function(j){
    Sigma$A[j] / Sigma$B[j]  * V_Lambda[,, j] %*% Eta$M %*% 
      (params$Z$M[, j] - X %*% params$Beta$M[, j, drop = FALSE]) 
  })
  list(M = matrix(M_Lambda, nrow = params$q, ncol = params$p), Cov = V_Lambda)
}


get_update_Poisson_VI_Eta <- function(params, X){
  # Useful quantities
  Lambda <- params$Lambda
  Sigma <- params$Sigma
  Phi <- params$Phi
  Delta <- params$Delta
  
  #First, get the common covariance matrix
  common_cov <- lapply(1:params$p, function(j){
    Sigma$A[j] / Sigma$B[j] * 
      (Lambda$M[,j] %*% t(Lambda$M[,j]) + Lambda$Cov[,,j])
  }) %>% 
    Reduce(f = "+") %>% 
    {. + diag(params$q)} %>% 
    solve()
  # The common matrix for all means
  mean_matrix <- common_cov %*% Lambda$M %*% diag(x = Sigma$A/Sigma$B, 
                                                  nrow = params$p, 
                                                  ncol = params$p)
  # The means
  # Response matrix
  response_matrix <- params$Z$M - X %*% params$Beta$M
  M <- sapply(1:params$n, function(i){
    mean_matrix %*% response_matrix[i, ] 
  })
  Cov <- array(common_cov, dim = c(params$q, params$q, params$n))
  list(M = matrix(M, nrow = params$q, ncol = params$n), Cov = Cov)
}

get_update_Poisson_VI_Sigma_without_fixed_effects <- function(params, priors){
  Lambda <- params$Lambda
  Eta <- params$Eta
  Phi <- params$Phi
  Delta <- params$Delta
  A <- rep(priors$Sigma$A + params$n / 2, params$p)
  E_eta_prime_eta <- apply(Eta$Cov, c(1, 2), sum) + # Sum of variances
    Reduce(f = "+", 
           x = lapply(1:params$n, function(i){
             Eta$M[, i] %*% t(Eta$M[, i])
           }))
  get_B_sigma_j <- function(j){
    term1 <- 0.5 * sum(params$Z$M[, j]^2  +  params$Z$S2[,j]) 
    term2 <- -sum(params$Z$M[, j] * (t(Eta$M) %*% Lambda$M[,j])) # Eta$M is coded in q x n
    term3 <- 0.5 * sum(E_eta_prime_eta * 
                         (Lambda$Cov[,, j] + Lambda$M[, j] %*% t(Lambda$M[, j]))) 
    term1 + term2 + term3
  }
  B <- priors$Sigma$B + sapply(1:params$p, get_B_sigma_j)
  list(A = A, B = B)
}

get_update_Poisson_VI_Sigma <- function(params, priors,
                                        X = 0, XprimeX = 0){
  Lambda <- params$Lambda
  Eta <- params$Eta
  Beta <- params$Beta
  # Posterior variationel de A
  update_without_X <- get_update_Poisson_VI_Sigma_without_fixed_effects(params, priors)
  A = update_without_X$A
  if(all(X == 0)){ # Case when no covariates
    return(list(A = A, B = update_without_X$B))
  }
  else{
    get_update_term_B_sigma_j <- function(j){
      term1 <- 0.5 * sum(XprimeX * (Beta$Cov[,,j] + Beta$M[,j] %*% t(Beta$M[,j])))
      term2 <- -sum((params$Z$M[, j] -  t(Eta$M) %*% Lambda$M[,j]) * (X %*% Beta$M[,j])) # Eta$M is coded in q x n
      term1 + term2
    }
    B <- update_without_X$B + map_dbl(1:params$p, get_update_term_B_sigma_j) 
  }
  return(list(A = A, B = B))
}

get_update_Poisson_VI_Z <- function(Y, X, params){
  target_function <- function(z_pars_ij, y_ij, A_j, B_j, M_ij) {
    # i de 1 à n, j de 1 à p
    m_ij = z_pars_ij[1]
    s2_ij = z_pars_ij[2]
    res = y_ij * m_ij - exp(m_ij + 0.5 * s2_ij)  - # Likelihood term 
      0.5 * A_j / B_j * (m_ij^2 + s2_ij - 2 * m_ij * M_ij) + 
      0.5 * log(s2_ij) # Variational entropy
    return(-res) # Return the negative results for optimization
  }
  # Matrix of prior expectations for Z
  prior_means <- X %*% params$Beta$M + t(params$Eta$M) %*% params$Lambda$M
  A <- params$Sigma$A
  B <- params$Sigma$B
  m_start <- params$Z$M
  s2_start <- params$Z$S2
  colnames(m_start) = NULL
  index_matrix <- expand.grid(i = 1:params$n,
                              j = 1:params$p)
  optim_results <- parallel::mcmapply(index_matrix[, "i"], 
                                      index_matrix[, "j"],
                                      FUN = function(i, j){
                                        Opt <- optim(
                                          par = c(m_start[i, j], s2_start[i, j]),
                                          fn = target_function,
                                          y_ij = Y[i, j],
                                          A_j = A[j],
                                          B_j = B[j],
                                          M_ij = prior_means[i,j],
                                          method = "L-BFGS-B",
                                          lower = c(-Inf, 1e-10),
                                          upper = c(Inf, 100)
                                        )
                                        c(mean = Opt$par[1], var = Opt$par[2])
                                      },
                                      mc.cores = parallel::detectCores() - 2,
                                      SIMPLIFY = TRUE)
  list(M = matrix(optim_results["mean", ],
                  nrow = params$n,
                  ncol = params$p),
       S2 = matrix(optim_results["var", ],
                   nrow = params$n,
                   ncol = params$p))
}

get_update_Poisson_amortized_VI_Z <- function(X, params){
  target_function <- function(Y, encoder, X, M, A, B) {
    # i de 1 à n, j de 1 à p
    res = 0
    for(i in 1:params$n){
    m_i = encoder(X)$M[i]
    s2_i = encoder(X)$S2[i]
    y_i = Y[i]
    M_i = M[i]
    res = res + y_i * m_i - exp(m_i + 0.5 * s2_i)  - # Likelihood term 
      0.5 * A / B * (m_i^2 + s2_i - 2 * m_i * M_i) + 
      0.5 * log(s2_i)
    }# Variational entropy
    return(torch_sum(-res)) # Return the negative results for optimization
  }
  # Matrix of prior expectations for Z
  prior_means <- X %*% params$Beta$M + t(params$Eta$M) %*% params$Lambda$M
  A <- params$Sigma$A
  B <- params$Sigma$B
  m_start <- params$encoder(X)$M
  s2_start <- params$encoder(X)$S2
  colnames(m_start) = NULL
  target = lapply(1:params$n, target_function)
  index_matrix <- expand.grid(i = 1:params$n,
                              j = 1:params$p)
  optim_results <- parallel::mcmapply(index_matrix[, "i"], 
                                      index_matrix[, "j"],
                                      FUN = function(i, j){
                                        Opt <- optim(
                                          par = c(m_start[i, j], s2_start[i, j]),
                                          fn = target_function,
                                          y_ij = Y[i, j],
                                          A_j = A[j],
                                          B_j = B[j],
                                          M_ij = prior_means[i,j],
                                          method = "L-BFGS-B",
                                          lower = c(-Inf, 1e-10),
                                          upper = c(Inf, 100)
                                        )
                                        c(mean = Opt$par[1], var = Opt$par[2])
                                      },
                                      mc.cores = parallel::detectCores() - 2,
                                      SIMPLIFY = TRUE)
  list(M = matrix(optim_results["mean", ],
                  nrow = params$n,
                  ncol = params$p),
       S2 = matrix(optim_results["var", ],
                   nrow = params$n,
                   ncol = params$p))
}



get_update_Poisson_VI_Beta <- function(params, priors, X, XprimeX){
  Sigma <- params$Sigma
  # Calcul des variances
  V_Beta <- lapply(1:params$p, function(j){
    precision <- Sigma$A[j] / Sigma$B[j] * XprimeX + priors$Beta$Precision[,,j]
    variance <- solve(precision)
    return(variance)
  }) %>% 
    abind(along = 3) # Mise au format array
  M_Beta <- sapply(1:params$p, function(j){
    V_Beta[,, j] %*% (
      Sigma$A[j] / Sigma$B[j]  * t(X) %*% 
        (params$Z$M[, j] - t(params$Eta$M) %*% params$Lambda$M[, j]) + 
        priors$Beta$Precision[,, j] %*% priors$Beta$M[,j])
  })
  list(M = matrix(M_Beta, nrow = params$F_x, ncol = params$p),
       Cov = V_Beta)
}

get_update_Poisson_VI_Phi <- function(params, priors){
  Lambda <- params$Lambda
  Delta <- params$Delta
  
  # le meme pour tous
  A <- matrix(priors$Phi$A + 0.5, params$p, params$q)
  
  cumprod_Delta  = cumprod(Delta$A/Delta$B)
  B <- sapply(1:params$q, function(h){
    priors$Phi$B + 
      .5 * (Lambda$M[h, ]^2 + Lambda$Cov[h, h, ]) * cumprod_Delta[h]
  })
  list(A = A, B = matrix(B, params$p, params$q))
}

get_update_Poisson_VI_Delta <- function(params, priors){
  new_A <- priors$Delta$A + 0.5 * params$p * (params$q + 1 - (1:params$q))
  E_phi_L2 <- sapply(1:params$q, function(h){
    E_phi_h <- params$Phi$A[, h] / params$Phi$B[, h]
    E_L2 <- params$Lambda$M[h, ]^2 + params$Lambda$Cov[h, h, ]
    sum(E_phi_h * E_L2)
  })
  new_B <- params$Delta$B
  new_B[1] <- priors$Delta$B + .5 * sum(E_phi_L2 * cumprod(new_A / new_B) / (new_A[1] / new_B[1]))
  if(params$q > 1){
    for(k in 2:params$q){
      E_delta_all <- cumprod(new_A / new_B)
      E_delta_k <- (new_A[k] / new_B[k])
      new_B[k] <- priors$Delta$B  + .5 * sum(E_phi_L2[k:params$q] * 
                                               E_delta_all[k:params$q] / E_delta_k )
    }
  }
  list(A = new_A, B = new_B)
}


# Miscellaneous functions for ELBO ----------------------------------------


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


# Main function to compute ELBO -------------------------------------------


get_ELBO <- function(Y, params, priors, X = 0, XprimeX = 0){
  # Variational entropy term
  variational_entropy <- sum(apply(params$Lambda$Cov, 3, 
                                   get_entropy_normal)) + # Lambda
    sum(apply(params$Eta$Cov, 3, get_entropy_normal)) + # Eta
    sum(.5 * log(params$Z$S2)) + # Somme des entropies pour Z
    sum(get_entropy_gamma(params$Sigma$A, params$Sigma$B)) + # Sigma 
    sum(get_entropy_gamma(params$Delta$A, params$Delta$B)) + # Delta 
    sum(get_entropy_gamma(params$Phi$A, params$Phi$B)) # Phi
  if(any(X != 0)){ # Case with covariates X with parameter beta
    variational_entropy <- variational_entropy +
      sum(apply(params$Beta$Cov, 3, get_entropy_normal))  # Beta
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
           x = lapply(1:params$n, function(i){
             params$Eta$M[, i] %*% t(params$Eta$M[, i])
           }))
  get_E_quadr_form <- function(j){
    term1 <- sum(0.5 * sum(params$Z$M[, j]^2 + params$Z$S2[, j]))
    term2 <- - sum(params$Z$M[, j] * (t(params$Eta$M) %*% params$Lambda$M[,j])) # Eta$M is coded in q x n
    term3 <- 0.5 * sum(E_eta_prime_eta * 
                         (params$Lambda$Cov[,, j] + params$Lambda$M[, j] %*% t(params$Lambda$M[, j]))) 
    term2bis <- 0
    term3bis <- 0
    term4 <- 0
    if(any(X != 0)){ # Case with covariates X with parameter beta
      term2bis <- -sum(params$Z$M[, j] * (X %*% params$Beta$M[,j]) ) 
      term3bis <- 0.5 * sum(XprimeX * 
                              (params$Beta$Cov[,, j] + params$Beta$M[, j] %*% t(params$Beta$M[, j]))) 
      term4 <- sum((X %*% params$Beta$M[,j])*t(params$Lambda$M[,j]%*%params$Eta$M))
    }
    
    term1 + term2 + term3 + term2bis + term3bis + term4
  }
  z_knowing_theta_expectation <- sum(.5 * params$n * expectations_log_sigma -
                                       expectations_sigma * map_dbl(1:params$p, get_E_quadr_form))
  poisson_log_likelihood <- sum(Y * params$Z$M - exp(params$Z$M + 0.5 * params$Z$S2))
  # Priors terms
  prior_sigma_expectation <- sum((priors$Sigma$A - 1) * expectations_log_sigma - 
                                   priors$Sigma$B * expectations_sigma)
  prior_phi_expectation <- sum((priors$Phi$A - 1) * expectations_log_phi - 
                                 priors$Phi$B * expectations_phi)
  prior_delta_expectation <- sum((priors$Delta$A - 1) * expectations_log_delta - 
                                   priors$Delta$B * expectations_delta)
  prior_eta_expectation <- -0.5 * sum(map_dbl(1:params$n, function(i){
    sum(params$Eta$M[,i]^2 + diag(x = matrix(params$Eta$Cov[,, i],
                                             nrow = params$q, ncol = params$q)))
  }))
  prior_lambda_expectation <- 0.5 * params$p * sum(cumsum(expectations_log_delta)) +
    0.5 * sum(expectations_log_phi) - 
    0.5 * sum(map_dbl(1:params$p, function(j){
      sum(cumprod(expectations_delta) * diag(x = expectations_phi[j, ], nrow = params$q, ncol = params$q) * 
            (diag(x = matrix(params$Lambda$Cov[,, j],
                             nrow = params$q, ncol = params$q)) + params$Lambda$M[, j]^2))
    }))
  prior_beta_expectation <- 0
  if(any(X != 0)){ 
    # On considère un prior normal (non lié aux sigmas)
    prior_beta_expectation <- 
      -0.5 * sum(map_dbl(1:params$p, function(j){
        # Expectation of beta'P0 beta
        sum(priors$Beta$Precision[,,j] * 
              (params$Beta$Cov[,, j] + diag(x = params$Beta$M[, j]^2, nrow = params$F_x, ncol = params$F_x))) -
          2 * sum(params$Beta$M[, j] * priors$Beta$Precision[,,j] %*% priors$Beta$M[,j])
        
      }))
  }
  ELBO <- variational_entropy +
    z_knowing_theta_expectation +
    poisson_log_likelihood +
    prior_delta_expectation +
    prior_eta_expectation +
    prior_lambda_expectation +
    prior_sigma_expectation +
    prior_phi_expectation +
    prior_beta_expectation
  return(ELBO)
}




get_CAVI <- function(Y, 
                     X = NULL, 
                     q, 
                     priors = NULL,
                     n_steps, 
                     seed = NULL,
                     params = NULL,
                     updates = NULL,
                     debug = FALSE,
                     get_ELBO_freq = 1){
  p <- ncol(Y); n <- nrow(Y); 
  # Checking priors
  if(is.null(priors)){
    print("A default prior was set, this should be avoided")
    priors = list(Sigma = list(A = 3, B = 2), 
                  Phi = list(A = 3/2, B = 3/2),
                  Delta= list(A = rep(3, q), 
                              B = 1))
    if(is.null(X)){
      priors$Beta = list(M = rep(0, 1),
                         C = rep(0.01, 1))
    }
    else{
      priors$Beta = list(M = matrix(0,
                                    nrow = ncol(X), ncol = p),
                         Precision = array(diag(0.01, ncol(X)),
                                           dim = c(ncol(X), ncol(X), p)))
    }
  }
  if(is.null(updates)){
    print("As no updates were provided, by default, all unknown are updated")
    updates = c(Lambda = TRUE, Sigma = TRUE,
                Eta = TRUE, Delta = TRUE, Phi = TRUE, 
                Beta = TRUE, Z = TRUE)
  }
  if(is.null(params)){
    print("As no initial parameters were fixed, there values is set randomly")
    if(!is.null(seed)){
      print(paste("Random seed fixed to", seed))
      set.seed(seed)
    }
    params <- list(Lambda = list(M = matrix(rnorm(p * q),
                                            nrow = q, ncol = p),
                                 Cov = array(diag(1, q), 
                                             dim = c(q, q, p))),
                   Eta = list(M = matrix(rnorm(n * q), 
                                         nrow = q, ncol = n),
                              Cov = array(diag(1, q), 
                                          dim = c(q, q, n))),
                   Sigma = list(A = rep(priors$Sigma$A + n / 2, p),
                                B = runif(p, 1, 5)),
                   Delta = list(A = priors$Delta$A + 0.5 * p * (q + 1 - (1:q)),
                                B = runif(q, 1, 10)),
                   Phi = list(A =  matrix(priors$Phi$A + 0.5, p, q),
                              B = matrix(runif(p * q, 1, 3), p, q)),
                   Z = list(M = log(Y + 1),
                            S2 = matrix(.1, nrow = nrow(Y),
                                        ncol = ncol(Y)))) 
    if(is.null(X)){
      params$Beta = list(M = matrix(0,
                                    nrow = 1, ncol = p),
                         Cov = array(diag(x = 0.001, nrow = 1),
                                     dim = c(1, 1, p)))
    }
    else{
      params$Beta = list(M = matrix(rnorm(p * ncol(X)),
                                    nrow = ncol(X), ncol = p),
                         Cov = array(diag(1, ncol(X)), 
                                     dim = c(ncol(X), ncol(X), p)))
    }
  }
  if(is.null(X)){
    updates["Beta"] <- FALSE
    print("As no X is provided, Beta won't be updated")
    X <- matrix(0, nrow = n, ncol = 1)
  }
  # Defining reccurrent objects
  XprimeX <- t(X) %*% X
  # Adding constant to params
  params$n <- nrow(Y)
  params$p <- ncol(Y)
  params$q <- q
  params$F_x <- ncol(X)
  
  current_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                           X = X, XprimeX = XprimeX)
  ELBOS <- data.frame(iteration = 0, 
                      ELBO = current_ELBO)
  pacman::p_load(progress)
  my_progress_bar <- progress_bar$new(total=n_steps)
  options(width = 80)
  for(step_ in 1:n_steps){
    
      my_progress_bar$tick()
    # if((step_ %% round(n_steps/10))==0){
    #   stepbar <- round(step_ / n_steps * (width - extra))
    #   text <- sprintf('|%s%s|% 3s%%', strrep('=', stepbar),
    #                   strrep(' ', width - stepbar - extra), round(step_ / n_steps * 100))
    #   cat(text)
    #   cat('\n')
    #   # Sys.sleep(0.05)
    #   # cat(if (step_ == n_steps) '\n' else '\014')
    #   #   my_progress_bar$tick()
    # }

    # Lambdas
    if(updates["Lambda"]){
      params$Lambda <- get_update_Poisson_VI_Lambda(params = params, X = X)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(new_ELBO - current_ELBO)
          print(paste("Problem at iteration", step_, "after updating Lambda"))
        }
        current_ELBO <- new_ELBO
      }
    }
    # Etas
    if(updates["Eta"]){
      params$Eta <- get_update_Poisson_VI_Eta(params = params, X = X)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(new_ELBO - current_ELBO)
          print(paste("Problem at iteration", step_, "after updating Eta"))
        }
        current_ELBO <- new_ELBO
      }
    }
    # Sigma
    if(updates["Sigma"]){
      params$Sigma <- get_update_Poisson_VI_Sigma(params = params, priors = priors,
                                                  X = X, XprimeX = XprimeX)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(new_ELBO - current_ELBO)
          print(paste("Problem at iteration", step_, "after updating Sigma"))
        }
        current_ELBO <- new_ELBO
      }
    }
    if(updates["Beta"]){
      params$Beta <- get_update_Poisson_VI_Beta(params = params, priors = priors, 
                                                X = X, XprimeX = XprimeX)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(new_ELBO - current_ELBO)
          print(paste("Problem at iteration", step_, "after updating Beta"))
        }
        current_ELBO <- new_ELBO
      }
    }
    # Phis
    if(updates["Phi"]){
      params$Phi <- get_update_Poisson_VI_Phi(params, priors)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print("Problem at iteration", step_, "after updating Phi")
        }
        current_ELBO <- new_ELBO
      }
    }
    if(updates["Delta"]){
      # deltas
      params$Delta <- get_update_Poisson_VI_Delta(params, priors)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(paste("Problem at iteration", step_, "after updating Delta"))
        }
        current_ELBO <- new_ELBO
      }
    }
    if(updates["Z"]){
      params$Z <- get_update_Poisson_VI_Z(Y = Y, X = X, params = params)
      if(debug){
        new_ELBO <- get_ELBO(Y = Y, params = params, priors = priors, 
                             X = X, XprimeX = XprimeX)
        if(new_ELBO < current_ELBO){
          print(new_ELBO - current_ELBO)
          print(paste("Problem at iteration", step_, "after updating Z"))
        }
        current_ELBO <- new_ELBO
      }
    }
    if((n_steps %% get_ELBO_freq) == 0){
      ELBOS <- bind_rows(ELBOS,
                         data.frame(iteration = step_,
                                    ELBO = get_ELBO(Y = Y, params = params, 
                                                    priors = priors, 
                                                    X = X, XprimeX = XprimeX)))
    }
  }
  return(list(ELBOS = ELBOS, params = params))
}



# Formatting functions ----------------------------------------------------



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
              pivot_longer(-c("true_paramsiteration", "row"), 
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