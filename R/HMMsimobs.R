sim_data_from_HMM <- function(N, T, init_state_probs, transition_matrix, obs_dist_params, robs_dist) {
  S <- length(init_state_probs)
  
  obs <- lapply(seq_len(N), function(i) {
    tempres <- sim_data_from_HMM_one_subject(T[i], init_state_probs, transition_matrix, obs_dist_params, robs_dist)
	tempres$subject <- i
	return(tempres)
  })

  return(obs)
}

sim_data_from_HMM_one_subject <- function(T, init_state_probs, transition_matrix, obs_dist_params, robs_dist) {
  state_space <- seq_along(init_state_probs)
  
  y <- rep(NA, T)
  
  y[1] <- sample(state_space, 1, prob = init_state_probs)
  for(t in seq_len(T - 1))
    y[t + 1] <- sample(state_space, 1, prob = transition_matrix[y[t],])
  
  X <- matrix(NA, nrow = T, ncol = obs_dist_params$D)
#  X <- matrix(NA, nrow = T, ncol = 8)
  
  for(s in state_space)
    X[y == s, ] <- robs_dist(sum(y == s), obs_dist_params[[s]])
  
  return(list(y = factor(y, levels = state_space), X = X))
}

sim_data_from_2HMM <- function(N, T, init_state_probs, transition_matrix, obs_dist_params, robs_dist) {
  S <- length(init_state_probs)
  
  obs <- lapply(seq_len(N), function(i) {
    tempres <- sim_data_from_2HMM_one_subject(T[i], init_state_probs, transition_matrix, obs_dist_params, robs_dist)
	tempres$subject <- i
	return(tempres)
  })

  return(obs)
}

sim_data_from_2HMM_one_subject <- function(T, init_state_probs, transition_matrix, obs_dist_params, robs_dist) {
  state_space <- seq_along(init_state_probs)
  
  y <- rep(NA, T)
  
  y[1] <- sample(state_space, 1, prob = init_state_probs)
  y[2] <- sample(state_space, 1, prob = apply(transition_matrix[, y[1], ], 2, sum))
  for(t in seq_len(T - 2))
    y[t + 2] <- sample(state_space, 1, prob = transition_matrix[y[t], y[t + 1],])
  
  X <- matrix(NA, nrow = T, ncol = obs_dist_params$D)
#  X <- matrix(NA, nrow = T, ncol = 8)
  
  for(s in state_space)
    X[y == s, ] <- robs_dist(sum(y == s), obs_dist_params[[s]])
  
  return(list(y = factor(y, levels = state_space), X = X))
}


sim_data_from_TDGMM <- function(N, T, init_state_probs, transition_matrix, duration_dist_params, obs_dist_params, robs_dist) {
  S <- length(init_state_probs)
  
  obs <- lapply(seq_len(N), function(i) {
    tempres <- sim_data_from_TDGMM_one_subject(T[i], init_state_probs, transition_matrix, duration_dist_params, obs_dist_params, robs_dist)
	tempres$subject <- i
	return(tempres)
  })

  return(obs)
}

sim_data_from_TDGMM_one_subject <- function(T, init_state_probs, transition_matrix, duration_dist_params, obs_dist_params, robs_dist) {
  state_space <- seq_along(init_state_probs)

  dd_start_state <- sapply(duration_dist_params, function(params_one_comp) params_one_comp$start_state)
  dd_end_state <- sapply(duration_dist_params, function(params_one_comp) params_one_comp$end_state)
  dd_head_probs <- lapply(duration_dist_params, function(params_one_comp) calc_duration_dist_head_probs(params_one_comp$alpha, params_one_comp$beta, params_one_comp$r, params_one_comp$M))
  
  y <- rep(NA, T)
  
  t <- 1
  y[t] <- sample(state_space, 1, prob = init_state_probs)
  while(t < T) {
    # sample next state given current state, based on transition matrix
    y[t + 1] <- sample(state_space, 1, prob = transition_matrix[y[t],])

	dd_ind <- which(dd_start_state == y[t] & dd_end_state == y[t + 1])

	# sample duration in next state, based on combination of states at times t and t + 1
	head_portion <- rbinom(1, size = 1, prob = duration_dist_params[[dd_ind]]$q)
	if(head_portion) {
		duration <- sample(seq_len(duration_dist_params[[dd_ind]]$M), size = 1, prob = dd_head_probs[[dd_ind]])
	} else {
		duration <- duration_dist_params[[dd_ind]]$M + 1 + rgeom(1, prob = 1 - duration_dist_params[[dd_ind]]$s)
	}

	# update y to the given state throughout its duration
	t <- t + 1
	t_offset <- 1
	while(t_offset < duration && t < T) {
		y[t + 1] <- y[t]
		t <- t + 1
		t_offset <- t_offset + 1
	}
  }
  
  X <- matrix(NA, nrow = T, ncol = obs_dist_params$D)
  X <- matrix(NA, nrow = T, ncol = 8)
  
  for(s in state_space)
    X[y == s, ] <- robs_dist(sum(y == s), obs_dist_params[[s]])
  
  return(list(y = factor(y, levels = state_space), X = X))
}


calc_duration_dist_head_probs <- function(alpha, beta, r, M) {
	probs <- sapply(seq_len(M), function(tau) {
		if(alpha == 0) {
			return(lgamma(alpha + beta) + lgamma(alpha + r) + lgamma(tau + r - 1) + lgamma(tau + beta - 1) - lgamma(r) - lgamma(beta) - lgamma(tau) - lgamma(tau + r + alpha + beta - 1))
		} else {
			return(lgamma(alpha + beta) + lgamma(alpha + r) + lgamma(tau + r - 1) + lgamma(tau + beta - 1) - lgamma(r) - lgamma(alpha) - lgamma(beta) - lgamma(tau) - lgamma(tau + r + alpha + beta - 1))
		}
	})

	probs <- exp(probs)

	return(probs / sum(probs))
}
