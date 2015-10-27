generate_transition_matrix_full_parameterization <- function(S = 2, alpha, log = TRUE) {
  if(missing(alpha))
    alpha <- lapply(seq_len(S), function(s) {
      temp <- rep(1, S)
      temp[s] <- 9
      return(temp)
    })
  
  trans_matrix <- matrix(NA, nrow = S, ncol = S)
  
  for(s in seq_len(S))
    trans_matrix[s, ] <- rdirichlet(n = 1, alpha = alpha[[s]])
  
  if(log) {
    return(log(trans_matrix))
  } else {
    return(trans_matrix)
  }
}

generate_transition_matrix_reduced_parameterization <- function(S = 2, alpha = 10, beta = 1.5, omega, log = TRUE) {
#  if(all(c(missing(alpha), missing(beta), missing(omega)))) {
#    stop("Either (both of alpha and beta) or omega must be supplied as arguments")
#  } else if(all(c(!missing(alpha), !missing(beta), !missing(omega)))) {
#    stop("All three of alpha, beta, and omega have been supplied.  Please supply either (alpha and beta) or omega, but not all three.")
#  }

  if(missing(omega)) {
    omega <- log(rbeta(1, alpha, beta))
  }
  
  if(omega > 0)
    stop("omega must be between -Inf and 0")
  
  log_trans_matrix <- matrix(logspace_sub(0, omega) - log(S - 1), nrow = S, ncol = S)
  diag(log_trans_matrix) <- omega

  if(log) {
    return(log_trans_matrix)
  } else {
    return(exp(log_trans_matrix))
  }
}

r_overlapping_GMM_params <- function(S, D, idim, M, alpha = rep(1, M), mu_scale = 5, Xi = rep(0, D), Psi = diag(D), Sigma_scale = 0.25, nu = D, iwish_scale = diag(D), cov_diag_expon = 0.5,
                                     rho_adj_var = 0.1 / D, shift_cov_scale = 0.25, cov_rotate_bounds, cov_rotate_var, cov_scale_var, cov_expon_bounds, cov_expon_mean, cov_expon_var) {
  rho <- rdirichlet(1, alpha)
  
  mus <- vector("list", M)
  Sigmas <- vector("list", M)
  for(m in seq_len(M)) {
    mus[[m]] <- as.numeric(mu_scale * rmvnorm(1, Xi, Psi))
    Sigmas[[m]] <- Sigma_scale * riwish(nu, iwish_scale)
    
    eSigma <- eigen(Sigmas[[m]])
    evals <- eSigma$values
    evecs <- eSigma$vectors
    
    diag_Sigma <- diag(evals^cov_diag_expon, nrow = D)
    
    Sigmas[[m]] <- evecs %*% diag_Sigma %*% t(evecs)
  }
  
  results <- vector("list", S)
  results[[1]] <- list(rho = rho, mus = mus, Sigmas = Sigmas)
  
  if(identical(idim, "F")) {
    translateD <- D
  } else if(identical(idim, "S")) {
    translateD <- 5L
  }
  
  mult_vec <- matrix(0, nrow = D, ncol = 1)
  for(s in 1 + seq_len(S - 1)) {
    new_rho <- rho
#    new_rho <- rho + rtmvnorm(1, mean = rep(0, M), sigma = rho_adj_var * diag(M), lower = -as.numeric(rho), upper = rep(1, M))
#    new_rho <- new_rho / sum(new_rho)
    
    new_mus <- mus
    new_Sigmas <- Sigmas
    
    for(m in seq_len(M)) {
      eSigma <- eigen(Sigmas[[m]])
      evals <- eSigma$values
      evecs <- eSigma$vectors
      
      final_rot_mat <- diag(translateD)
      for(i in seq_len(translateD - 1)) {
        for(j in seq(from = i + 1, to = translateD, by = 1)) {
          theta <- rtmvnorm(1, mean = 0, sigma = matrix(cov_rotate_var), cov_rotate_bounds$lower, cov_rotate_bounds$upper)
          cos_theta <- cos(theta)
          sin_theta <- sin(theta)
          
          Givens_rot_mat <- diag(translateD)
          Givens_rot_mat[i, i] <- cos_theta
          Givens_rot_mat[j, j] <- cos_theta
          Givens_rot_mat[i, j] <- -1 * sin_theta
          Givens_rot_mat[j, i] <- sin_theta
          
          final_rot_mat <- final_rot_mat %*% Givens_rot_mat
        }
      }
      
      new_evecs <- evecs
      new_evecs[seq_len(translateD), seq_len(translateD)] <- evecs[seq_len(translateD), seq_len(translateD)] %*% final_rot_mat
      new_cov_scale <- as.numeric(rtmvnorm(1, mean = 1, sigma = matrix(cov_scale_var), lower = 0))
      new_cov_diag_expon <- as.numeric(rtmvnorm(1, mean = cov_expon_mean, sigma = matrix(cov_expon_var), lower = cov_expon_bounds$lower, upper = cov_expon_bounds$upper))
      new_Sigmas[[m]] <- new_cov_scale * new_evecs %*% diag(c(evals[seq_len(translateD)]^new_cov_diag_expon, evals[translateD + seq_len(D - translateD)]), nrow = D) %*% t(new_evecs)
      
      new_mus[[m]][1:translateD] <- rmvnorm(1, mean = mus[[m]][1:translateD], sigma = shift_cov_scale * new_Sigmas[[m]][1:translateD, 1:translateD, drop = FALSE])
    }
    
    results[[s]] <- list(rho = new_rho, mus = new_mus, Sigmas = new_Sigmas)
  }
  
  return(results)
}

