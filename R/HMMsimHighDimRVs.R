rHMMSimStudyObsDistFulldim <- function(N, obs_dist_params) {

}

rHMMSimStudyObsDistSmalldim <- function(N, obs_dist_params) {
  D <- ncol(obs_dist_params$Sigmas[[1]])
  
  reduced_mus <- lapply(obs_dist_params$mus, function(mu) mu[1:5])
  reduced_Sigmas <- lapply(obs_dist_params$Sigmas, function(Sigma) Sigma[1:5, 1:5])
  
  return(cbind(rGMM(N, rho, reduced_mus, reduced_Sigmas),
               matrix(rnorm(N * (D - 5), 0, 100^2), nrow = N, ncol = (D - 5))))
}

dSimStudyDataUnbddFulldim <- function(X, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, log) {
  D <- ncol(X)
  
  if(D < 3)
    stop("D must be >= 3")
  
  tree_as_tnorms <- tree_to_tnorms(tree, tree_rotation_matrix)
  
  reduced_mus <- lapply(mus, function(mu) mu[1:(D - 2)])
  reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:(D - 2), 1:(D - 2)])
  
  rot_X <- X %*% t(overall_rotation_matrix)
  
  log_res <- apply(cbind(dTGMM(rot_X[, 1:2]  %*% tree_as_tnorms$rotation_matrix,
                               tree_as_tnorms$rho, tree_as_tnorms$mus, tree_as_tnorms$Sigmas, tree_as_tnorms$lower, tree_as_tnorms$upper, log = FALSE),
                         dGMM(rot_X[, 3:D], rhos = rho, mus = reduced_mus, Sigmas = reduced_Sigmas, log = TRUE)),
                   1,
                   sum)
  
  if(log) {
    return(log_res)
  } else {
    return(exp(log_res))
  }
}

dSimStudyDataUnbddSmalldim <- function(X, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, log) {
  D <- ncol(X)
  
  if(D < 5)
    stop("D must be >= 5")
  
  tree_as_tnorms <- tree_to_tnorms(tree, tree_rotation_matrix)
  
  reduced_mus <- lapply(mus, function(mu) mu[1:3])
  reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:3, 1:3])
  
  rot_X <- X %*% t(overall_rotation_matrix)
  
  temp <- cbind(dTGMM(rot_X[, 1:2]  %*% tree_as_tnorms$rotation_matrix,
                      tree_as_tnorms$rho, tree_as_tnorms$mus, tree_as_tnorms$Sigmas, tree_as_tnorms$lower, tree_as_tnorms$upper, log = FALSE),
                dGMM(rot_X[, 3:5], rhos = rho, mus = reduced_mus, Sigmas = reduced_Sigmas, log = TRUE))
  
  if(D > 5) {
    temp <- cbind(temp, dmvnorm(rot_X[, 6:D, drop = FALSE], rep(0, D - 5), diag(100^2, nrow = D - 5), log = TRUE))
  }
  
  log_res <- apply(temp, 1, sum)
  
  if(log) {
    return(log_res)
  } else {
    return(exp(log_res))
  }
}

rOverlappingGMMParams <- function(S, D, M, alpha = rep(1, M), mu_scale = 5, Xi = rep(0, D), Psi = diag(D), Sigma_scale = 0.25, nu = D, iwish_scale = diag(D),
                                  rho_adj_var = 0.1 / D, shift_cov_scale = 0.25, cov_rotate_D = ceiling(D / 3), cov_rotate_bounds, cov_rotate_var, cov_scale_var) {
  rho <- rdirichlet(1, alpha)
  
  mus <- vector("list", M)
  Sigmas <- vector("list", M)
  for(m in seq_len(M)) {
    mus[[m]] <- as.numeric(mu_scale * rmvnorm(1, Xi, Psi))
    Sigmas[[m]] <- Sigma_scale * riwish(nu, iwish_scale)
  }
  
  results <- vector("list", S)
  results[[1]] <- list(rho = rho, mus = mus, Sigmas = Sigmas)
  
#  Givens_rot_mat <- diag(D)
  mult_vec <- matrix(0, nrow = D, ncol = 1)
  for(s in 1 + seq_len(S - 1)) {
    new_rho <- rho + rtmvnorm(D, mean = 0, sigma = matrix(rho_adj_var))
    new_rho <- new_rho / sum(new_rho)
    
    new_mus <- mus
    new_Sigmas <- Sigmas
    
    for(m in seq_len(M)) {
      new_mus[[m]] <- rmvnorm(1, mean = mus[[m]], sigma = shift_cov_scale * Sigmas[[m]])
      
      for(i in seq_len(D - 1)) {
        for(j in seq(from = i + 1, to = D, by = 1)) {
          theta <- rtmvnorm(1, mean = 0, sigma = matrix(cov_rotate_var), cov_rotate_bounds$lower, cov_rotate_bounds$upper)
          cos_theta <- cos(theta)
          sin_theta <- sin(theta)
          
          mult_vec[i, 1] <- cos_theta
          mult_vec[j, 1] <- sin_theta
          new_Sigmas[[m]][, i] <- new_Sigmas[[m]] %*% mult_vec
          
          mult_vec[i, 1] <- -1 * sin_theta
          mult_vec[j, 1] <- cos_theta
          new_Sigmas[[m]][, j] <- new_Sigmas[[m]] %*% mult_vec
          
          mult_vec[i, 1] <- 0
          mult_vec[j, 1] <- 0
          
#          Givens_rot_mat[i, i] <- cos_theta
#          Givens_rot_mat[j, j] <- cos_theta
#          Givens_rot_mat[i, j] <- -1 * sin_theta
#          Givens_rot_mat[j, i] <- sin_theta
        }
      }
      
      cov_scale <- rnorm(1, 1, sqrt(cov_scale_var))
      new_Sigmas[[m]] <- cov_scale * new_Sigmas[[m]]
    }
    
    results[[s]] <- list(rho = new_rho, mus = new_mus, Sigmas = new_Sigmas)
  }
  
  return(results)
}
