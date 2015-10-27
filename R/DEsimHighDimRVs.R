rSimStudyDataUnbddFulldim <- function(N, D, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas) {
	reduced_mus <- lapply(mus, function(mu) mu[1:(D - 2)])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:(D - 2), 1:(D - 2)])
  
	tree_as_tnorms <- tree_to_tnorms(tree, tree_rotation_matrix)

	return(cbind(rTGMM(N, tree_as_tnorms$rho, tree_as_tnorms$mus, tree_as_tnorms$Sigmas, tree_as_tnorms$lower, tree_as_tnorms$upper) %*% t(tree_as_tnorms$rotation_matrix),
		rGMM(N, rho, reduced_mus, reduced_Sigmas)) %*% overall_rotation_matrix)
}

rSimStudyDataUnbddSmalldim <- function(N, D, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas) {
	reduced_mus <- lapply(mus, function(mu) mu[1:3])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:3, 1:3])
	
	tree_as_tnorms <- tree_to_tnorms(tree, tree_rotation_matrix)
	
	return(cbind(rTGMM(N, tree_as_tnorms$rho, tree_as_tnorms$mus, tree_as_tnorms$Sigmas, tree_as_tnorms$lower, tree_as_tnorms$upper) %*% t(tree_as_tnorms$rotation_matrix),
		rGMM(N, rho, reduced_mus, reduced_Sigmas),
		matrix(rnorm(N * (D - 5), 0, 100), nrow = N, ncol = (D - 5))) %*% overall_rotation_matrix)
}

rSimStudyDataBddFulldim <- function(N, D, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, lower, upper) {
	reduced_mus <- lapply(mus, function(mu) mu[1:(D - 2)])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:(D - 2), 1:(D - 2)])
	reduced_lower <- lapply(lower, function(lb) lb[1:(D - 2)])
	reduced_upper <- lapply(upper, function(ub) ub[1:(D - 2)])
	
	return(cbind(rpdf_tree(N, tree_rotation_matrix, tree, max_depth = NULL),
					rTGMM(N, rho, reduced_mus, reduced_Sigmas, reduced_lower, reduced_upper)) %*% overall_rotation_matrix)
}

rSimStudyDataBddSmalldim <- function(N, D, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, lower, upper) {
	reduced_mus <- lapply(mus, function(mu) mu[1:3])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:3, 1:3])
	reduced_lower <- lapply(lower, function(lb) lb[1:3])
	reduced_upper <- lapply(upper, function(ub) ub[1:3])
	
	return(cbind(rpdf_tree(N, tree_rotation_matrix, tree, max_depth = NULL),
					rTGMM(N, rho, reduced_mus, reduced_Sigmas, reduced_lower, reduced_upper),
					matrix(runif(N * (D - 5), 0, 100), nrow = N, ncol = (D - 5))) %*% overall_rotation_matrix)
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

dSimStudyDataBddFulldim <- function(X, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, lower, upper, log) {
	D <- ncol(X)
	
	if(D < 3)
		stop("D must be >= 3")
	
	reduced_mus <- lapply(mus, function(mu) mu[1:(D - 2)])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:(D - 2), 1:(D - 2), drop = FALSE])
	
	rot_X <- X %*% t(overall_rotation_matrix)
	
	temp <- cbind(dpdf_tree(rot_X[, 1:2], rotation_matrix = tree_rotation_matrix, tree = tree, max_depth = NULL, log = TRUE),
		dTGMM(rot_X[, 3:D], rho, reduced_mus, reduced_Sigmas, lower, upper, log = TRUE))
	
	temp <- apply(temp, 1, sum)
	
	if(log) {
		return(temp)
	} else {
		return(exp(temp))
	}
}

dSimStudyDataBddSmalldim <- function(X, overall_rotation_matrix, tree_rotation_matrix, tree, rho, mus, Sigmas, lower, upper, log) {
	D <- ncol(X)
	
	if(D < 5)
		stop("D must be >= 5")
	
	reduced_mus <- lapply(mus, function(mu) mu[1:3])
	reduced_Sigmas <- lapply(Sigmas, function(Sigma) Sigma[1:3, 1:3])
	reduced_lower <- lapply(lower, function(lb) lb[1:3])
	reduced_upper <- lapply(upper, function(ub) ub[1:3])
	
	rot_X <- X %*% t(overall_rotation_matrix)
	
	temp <- cbind(dpdf_tree(rot_X[, 1:2], rotation_matrix = tree_rotation_matrix, tree = tree, max_depth = NULL, log = TRUE),
		dTGMM(rot_X[, 3:5], rho, reduced_mus, reduced_Sigmas, reduced_lower, reduced_upper, log = TRUE))
	
	if(D > 5) {
	  temp <- cbind(temp, matrix(dunif(as.numeric(rot_X[, 6:D, drop = FALSE]), 0, 100, log = TRUE), ncol = D - 5))
	}
  
	temp <- apply(temp, 1, sum)
	
	if(log) {
		return(temp)
	} else {
		return(exp(temp))
	}
}


rGammaMM <- function(N = 1, rho = 1, component_params = list(component1 = list(var1 = list(gamma_loc = c(3), gamma_alpha = c(1), gamma_beta = c(2))))) {
	# rho is a vector of length M, where M is the number of mixture components
	# component_params is a list of length M with parameters for each component
	# element m of component_params is a list of length D, where D is the number of covariates to generate
	# component_params[[m]][[d]] has the parameters used to generate the d'th covariate from the m'th mixture component:
	#  - gamma_loc, a vector of length d used to compute location offsets
	#  - gamma_alpha, a vector of length d used to compute alpha
	#  - gamma_beta, a vector of length d used to compute beta

	M <- length(rho)
	
	D <- length(component_params[[1]])
	m <- sample.int(M, size = N, replace = TRUE, prob = rho)
	X <- matrix(NA, nrow = N, ncol = D)

	for(i in seq_len(M)) {
		if(sum(m == i) > 0)
			X[m == i, ] <- rmvGamma(sum(m == i), component_params = component_params[[i]])
	}
	
	return(X)
}

rmvGamma <- function(N, component_params) {
	result <- matrix(seq_len(0), nrow = N)

	for(var_ind in seq_along(component_params)) {
		var_params <- component_params[[var_ind]]

		if(var_ind == 1) {
			loc <- rep(var_params$gamma_loc[1], N)
			alpha <- rep(var_params$gamma_alpha[1], N)
			beta <- rep(var_params$gamma_beta[1], N)
		} else {
			loc <- var_params$gamma_loc[1] + result %*% matrix(var_params$gamma_loc[1 + seq_len(length(var_params$gamma_loc) - 1)])
			alpha <- var_params$gamma_alpha[1] + result %*% matrix(var_params$gamma_alpha[1 + seq_len(length(var_params$gamma_alpha) - 1)])
			beta <- var_params$gamma_beta[1] + result %*% matrix(var_params$gamma_beta[1 + seq_len(length(var_params$gamma_beta) - 1)])
		}

		result <- cbind(result, loc +  mapply(function(alpha_i, beta_i) { rgamma(1, shape = alpha_i, scale = beta_i) }, alpha, beta))

		if(any(result[, ncol(result)] == loc)) {
			browser()
		}
	}

	return(result)
}

dGammaMM <- function(X, rho = 1, component_params = list(component1 = list(var1 = list(gamma_loc = c(3), gamma_alpha = c(1), gamma_beta = c(2)))), log) {
	# rho is a vector of length M, where M is the number of mixture components
	# component_params is a list of length M with parameters for each component
	# element m of component_params is a list of length D, where D is the number of covariates to generate
	# component_params[[m]][[d]] has the parameters used to generate the d'th covariate from the m'th mixture component:
	#  - gamma_loc, a vector of length d used to compute location offsets
	#  - gamma_alpha, a vector of length d used to compute alpha
	#  - gamma_beta, a vector of length d used to compute beta

	M <- length(rho)
	D <- length(component_params[[1]])

	temp <- sapply(seq_len(M), function(m) {
		return(dmvGamma(X, component_params = component_params[[m]], log = TRUE))
	})

	temp <- log(rho) + temp
	res <- logspace_sum_matrix_rows(matrix(temp, nrow = 1))

	if(log) {
		return(res)
	} else {
		return(exp(res))
	}
}

dmvGamma <- function(X, component_params, log = TRUE) {
	if(is.vector(X))
		X <- matrix(X, nrow = 1)

	N <- nrow(X)

	result <- 0

	for(var_ind in seq_along(component_params)) {
		var_params <- component_params[[var_ind]]

		if(var_ind == 1) {
			loc <- rep(var_params$gamma_loc[1], N)
			alpha <- rep(var_params$gamma_alpha[1], N)
			beta <- rep(var_params$gamma_beta[1], N)
		} else {
			loc <- var_params$gamma_loc[1] + X[, seq_len(var_ind - 1), drop = FALSE] %*% matrix(var_params$gamma_loc[1 + seq_len(length(var_params$gamma_loc) - 1)])
			alpha <- var_params$gamma_alpha[1] + X[, seq_len(var_ind - 1), drop = FALSE] %*% matrix(var_params$gamma_alpha[1 + seq_len(length(var_params$gamma_alpha) - 1)])
			beta <- var_params$gamma_beta[1] + X[, seq_len(var_ind - 1), drop = FALSE] %*% matrix(var_params$gamma_beta[1 + seq_len(length(var_params$gamma_beta) - 1)])
		}

		result <- result + sum(sapply(seq_len(nrow(X)), function(obs_ind) { dgamma(X[obs_ind, var_ind] - loc[obs_ind], shape = alpha[obs_ind], scale = beta[obs_ind], log = TRUE) }))
	}

	if(log) {
		return(result)
	} else {
		return(exp(result))
	}
}

rGMM <- function(N = 1, rho = 1, mus = list(0), Sigmas = list(matrix(1))) {
	M <- length(rho)
	
	D <- length(mus[[1]])
	m <- sample.int(M, size = N, replace = TRUE, prob = rho)
	X <- matrix(NA, nrow = N, ncol = D)

	for(i in seq_len(M)) {
		if(sum(m == i) > 0)
			X[m == i, ] <- rmvnorm(sum(m == i), mean = mus[[i]], sigma = Sigmas[[i]])
#	return(mvrnorm(1, mu = model_params[[ind]]$mu, Sigma = model_params[[ind]]$Sigma))
	}
	
	return(X)
}

rTGMM <- function(N, rho, mus, Sigmas, lower, upper, constraints) {
	M <- length(mus)
	
	D <- length(mus[[1]])
	m <- sample.int(M, size = N, replace = TRUE, prob = rho)
	X <- matrix(NA, nrow = N, ncol = D)

	for(i in seq_len(M)) {
		if(sum(m == i) > 0) {
			if(missing(constraints)) {
				X[m == i, ] <- rtmvnorm(sum(m == i), mean = mus[[i]], sigma = Sigmas[[i]], lower[[i]], upper[[i]])
			} else {
				X[m == i, ] <- rtmvnorm(sum(m == i), mean = mus[[i]], sigma = Sigmas[[i]], lower[[i]], upper[[i]], D = constraints[[i]])
			}
		}
#	return(mvrnorm(1, mu = model_params[[ind]]$mu, Sigma = model_params[[ind]]$Sigma))
	}
	
	return(X)
}

dGMMR <- function(x, params, log = FALSE) {
	x <- as.matrix(x)
	res <- apply(x, 1, function(x_one) {
		return(dGMMR_one_obs(x = x_one, M = params$M, rho = params$rho, model_params = params$model_params, log = log))
	})
	
	return(res)
}

dTGMM <- function(X, rho, mus, Sigmas, lower, upper, constraints = NULL, norm_consts = rep(NA, length(rho)), log) {
	M <- length(rho)
	
	if(M > 0) {
		res <- dTnorm(X, mean = mus[[1]], sigma = Sigmas[[1]],
			lower = lower[[1]], upper = upper[[1]], constraints = constraints[[1]], norm_const = NULL, log = TRUE) + log(rho[1])
#		dtmvt(X, mean = mus[[1]], sigma = Sigmas[[1]], df = 0, lower = lower[[1]], upper = upper[[1]], log = TRUE)
		
		for(m in (seq_len(M - 1) + 1)) {
			res <- logspace_sum_matrix_rows(cbind(res,
							dTnorm(X, mean = mus[[m]], sigma = Sigmas[[m]],
								lower = lower[[m]], upper = upper[[m]], constraints = constraints[[m]], norm_const = NULL, log = TRUE) + log(rho[m])))
#							dtmvt(X, mean = mus[[m]], sigma = Sigmas[[m]], df = 0, lower = lower[[m]], upper = upper[[m]], log = TRUE)))
		}
	}
	
	return(res)
}

dTnorm <- function(X, mean, sigma, lower, upper, constraints = NULL, norm_const = NULL, log) {
	if(is.vector(X))
		X <- matrix(X)

	if(!is.null(constraints)) {
        mean <- as.vector(constraints %*% matrix(mean))
        sigma <- constraints %*% sigma %*% t(constraints)
        X <- X %*% t(constraints)
	}
	
	if(is.null(norm_const))
		norm_const <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)
	
	temp <- apply(X, 1,
		function(row) {
			all(c(row >= lower, row <= upper))
		})
	
	res <- rep(-Inf, nrow(X))
  
	if(sum(temp) > 0) {
		res[temp] <- dmvnorm(X[temp, , drop = FALSE], mean, sigma, log = TRUE) - log(norm_const)

		if(!is.null(constraints)) {
			res[temp] <- res[temp] + determinant(constraints)$modulus
		}
	}
	
	if(log) {
		return(res)
	} else {
		return(exp(res))
	}
}

dGMMR_one_obs <- function(x, M, rho, model_params, log = FALSE) {
	temp <- sapply(seq_len(M), function(m) {
		return(dmvnorm(x, mean = model_params[[m]]$mu, sigma = model_params[[m]]$Sigma, log = TRUE))
	})

	temp <- log(rho) + temp
	res <- logspace_sum_matrix_rows(matrix(temp, nrow = 1))

	if(log) {
		return(res)
	} else {
		return(exp(res))
	}
}

rGMMParams <- function(D, M, alpha = rep(1, M), mu_scale = 5, Xi = rep(0, D), Psi = diag(D), Sigma_scale = 0.25, nu = D, S = diag(D)) {
	rho <- rdirichlet(1, alpha)
	
	mus <- vector("list", M)
	Sigmas <- vector("list", M)
	for(m in seq_len(M)) {
		mus[[m]] <- as.numeric(mu_scale * rmvnorm(1, Xi, Psi))
		Sigmas[[m]] <- Sigma_scale * riwish(nu, S)
	}
	
	return(list(rho = rho, mus = mus, Sigmas = Sigmas))
}

getTGMMParams <- function(GMMParams) {
	M <- length(GMMParams$mus)
	D <- length(GMMParams$mus[[1]])
	
	target_p_center <- 0.5
	
	lower <- c()
	upper <- c()
	for(m in seq_len(M)) {
		marginal_p_tail <- optimize(opt_marginal_p_tail_difference, interval = c((1 - target_p_center)/(2*D), 0.5), target_p_center,
			mean = GMMParams$mus[[m]], Sigma = GMMParams$Sigmas[[m]], tol = 10^-4)
		
		temp <- qnorm(1 - marginal_p_tail$minimum) * sqrt(diag(GMMParams$Sigmas[[m]]))
		lower <- cbind(lower, GMMParams$mus[[m]] - temp)
		upper <- cbind(upper, GMMParams$mus[[m]] + temp)
	}
	
	lower <- apply(lower, 1, min)
	upper <- apply(upper, 1, max)
	
	M <- length(GMMParams$rho)
	lower <- lapply(seq_len(M), function(m) lower)
	upper <- lapply(seq_len(M), function(m) upper)
	
	return(list(rho = GMMParams$rho, mus = GMMParams$mus, Sigmas = GMMParams$Sigmas, lower = lower, upper = upper))
}

opt_marginal_p_tail_difference <- function(marginal_p_tail, target_p_center, mean, Sigma) {
	temp <- qnorm(1 - marginal_p_tail) * sqrt(diag(Sigma))
	(pmvnorm(lower = mean - temp, upper = mean + temp, mean = mean, sigma = Sigma) - target_p_center)^2
}


tree_to_tnorms <- function(tree, rotation_matrix) {
	temp <- get_leaf_bounds_and_log_probs(tree)
	
	lower <- lapply(temp$bounds, function(bds) {
		temp <- bds[1, ]
		temp[temp == tree$root_lower_endpts] <- -Inf
		return(temp)
	})
	
	upper <- lapply(temp$bounds, function(bds) {
		temp <- bds[2, ]
		temp[temp == tree$root_upper_endpts] <- Inf
		return(temp)
	})
	
	rho <- exp(temp$log_probs)
	
	mus <- lapply(temp$bounds, function(bds) {
		apply(bds, 2, mean)
	})
	
	Sigmas <- lapply(temp$bounds, function(bds) {
    diag((bds[2, ] - bds[1, ])/6, nrow = ncol(bds))
	})
	
	return(list(rho = rho, mus = mus, Sigmas = Sigmas, lower = lower, upper = upper, rotation_matrix = rotation_matrix))
}

unif_orthogonal_matrix_subgroup <- function(n = 2L, rotation_target_dimension = NULL) {
  if(n < 2)
    stop("n must be >= 2")
  
  n <- as.integer(n)
  
  result <- matrix(NA, nrow = n, ncol = n)
  
  if(identical(as.integer(n), 2L)) {
    multiplier <- NA
    if(is.null(rotation_target_dimension)) {
      multiplier <- sample(c(-1, 1), size = 1)
    } else if(identical(rotation_target_dimension %% 2, 0)) {
      multiplier <- 1
    } else if(identical(rotation_target_dimension %% 2, 1)) {
      multiplier <- -1
    }
    
    theta <- runif(1, 0, 2*pi)
    costheta <- cos(theta)
    sintheta <- sin(theta)
    result <- matrix(c(costheta, multiplier * sintheta, -1 * sintheta, multiplier * costheta), nrow = 2L, ncol = 2L)
  } else {
    z <- rnorm(n)
    v <- z/sqrt(sum(z^2))
    e <- c(1, rep(0, n - 1))
    x <- e - v
    x <- x/sqrt(sum(x^2))
    result_nm1 <- unif_orthogonal_matrix_subgroup(n - 1L, rotation_target_dimension)
    temp <- rbind(c(1, rep(0, n - 1)), cbind(rep(0, n - 1), result_nm1))
    result <- (diag(n) - 2 * outer(x, x)) %*% temp
  }
  
  return(result)
}

unif_rotation_matrix <- function(n) {
  return(unif_orthogonal_matrix_subgroup(n, rotation_target_dimension = n))
}

estimate_KL_divergence <- function(sample_true_log_density, sample_est_log_density) {
	if(identical(length(sample_true_log_density), length(sample_est_log_density))) {
		return( (sum(sample_true_log_density) - sum(sample_est_log_density))/length(sample_true_log_density) )
	} else {
		stop("True log density and estimated log density don't have the same length.")
	}
}

estimate_Hellinger_dist <- function(sample_true_log_density, sample_est_log_density) {
  n <- length(sample_true_log_density)
	if(identical(n, length(sample_est_log_density))) {
#  	temp <- mean(sqrt(exp(sample_est_log_density - sample_true_log_density)))
#	  if(temp > 1) {
#	    return(temp)
#	  } else {
#	    return( sqrt(1 - temp) )
#	  }
	  temp <- logspace_sum_matrix_rows(matrix(0.5 * (sample_est_log_density - sample_true_log_density) - log(n), nrow = 1))
	  if(temp > 0) {
			return(temp)
		} else {
			return( 0.5 * logspace_sub(0, temp) )
		}
	} else {
		stop("True log density and estimated log density don't have the same length.")
	}
}


estimate_Hellinger_dist_symmetric <- function(true_sample_true_log_density, true_sample_est_log_density, est_sample_true_log_density, est_sample_est_log_density) {
	if(!identical(length(true_sample_true_log_density), length(true_sample_est_log_density)) || 
		!identical(length(est_sample_true_log_density), length(est_sample_est_log_density))) {
		stop("True log density and estimated log density don't have the same length.")
	} else {
		temp <- c(estimate_Hellinger_dist(true_sample_true_log_density, true_sample_est_log_density),
				estimate_Hellinger_dist(est_sample_est_log_density, est_sample_true_log_density))
		if(all(temp < 1)) {
			return(mean(temp))
		} else {
			return(temp[2])
		}
	}
}


estimate_H_dist_D <- function(D, sample_true_log_density, sample_est_log_density) {
	if(!identical(length(sample_true_log_density), length(sample_est_log_density))) {
		stop("True log density and estimated log density don't have the same length.")
	}
	
	if(identical(as.integer(D), 1L) || identical(as.integer(D %% 2), 0L)) {
		temp <- (1 / D) * (sample_est_log_density - sample_true_log_density)
		
		for(i in seq_along(temp)) {
			if(temp[i] >= 0) {
				temp[i] <- D * logspace_sub(temp[i], 0)
			} else {
				temp[i] <- D * logspace_sub(0, temp[i])
			}
		}
		
		return((logspace_sum_matrix_rows(matrix(temp, nrow = 1)) - log(length(sample_true_log_density)) - log(D)) / D)
#		return(((1 / D) * exp(logspace_sum_matrix_rows(matrix(temp, nrow = 1))) / length(sample_true_log_density))^(1 / D))
	} else {
		stop("D must be equal to 1 or even for this function to work")
	}
}

estimate_H_dist_D_symmetric <- function(D, true_sample_true_log_density, true_sample_est_log_density, est_sample_true_log_density, est_sample_est_log_density) {
	if(!identical(length(true_sample_true_log_density), length(true_sample_est_log_density)) || 
		!identical(length(est_sample_true_log_density), length(est_sample_est_log_density))) {
		stop("True log density and estimated log density don't have the same length.")
	}
	
	return(mean(c(estimate_H_dist_D(D, true_sample_true_log_density, true_sample_est_log_density),
		estimate_H_dist_D(D, est_sample_est_log_density, est_sample_true_log_density))))
}

estimate_ISE <- function(true_sample_true_log_density, true_sample_est_log_density, est_sample_true_log_density, est_sample_est_log_density) {
	n1 <- length(true_sample_true_log_density)
	if(!identical(n1, length(true_sample_est_log_density)))
		stop("true_sample_true_log_density and true_sample_est_log_density must have the same length")
	
	n2 <- length(est_sample_true_log_density)
	if(!identical(n2, length(est_sample_est_log_density)))
		stop("est_sample_true_log_density and est_sample_est_log_density must have the same length")
	
#	term1 <- logspace_sum_matrix_rows(rbind(true_sample_true_log_density, true_sample_est_log_density)) - log(n1)
#	term2 <- logspace_sum_matrix_rows(rbind(est_sample_true_log_density, est_sample_est_log_density)) - log(n2)
#	
#	a <- logspace_add(term1[1], term2[1])
#	b <- logspace_add(term1[2], term2[2])
#	
#	return(exp(logspace_sub(a, b)))

	term1 <- logspace_sum_matrix_rows(rbind(true_sample_true_log_density, true_sample_est_log_density)) - log(n1)
	term2 <- logspace_sum_matrix_rows(rbind(est_sample_true_log_density, est_sample_est_log_density)) - log(n2)
	
#	a <- logspace_add(term1[1], term2[1])
#	b <- logspace_add(term1[2], term2[2])
	a <- logspace_add(term1[1], term2[2])
	b <- logspace_add(term1[2], term2[1])
	
  return(logspace_sub(a, b))
#	return(exp(logspace_sub(a, b)))
}
