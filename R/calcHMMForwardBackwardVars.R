calc_HMM_forward_vars_given_log_obs_probs <- function(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = TRUE) {
  # S is the number of elements in the state space
  # initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
  # trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
  
  state_space <- matrix(seq_len(S))
  
  forward_vars <- vector("list", T)
  
  forward_vars[[1]] <- matrix(log_initial_state_distn + log_obs_probs[, 1], nrow = 1)
  for(t in 2:T) {
    forward_vars[[t]] <- forward_vars[[t-1]]
    
    for(s in 1:S) {
      temp <- forward_vars[[t-1]] + log_trans_matrix[, s] + log_obs_probs[s, t]
      dim(temp) <- c(1, S)
      
      forward_vars[[t]][1, s] <- logspace_sum_matrix_rows(temp)
    }
  }
  
  if(!log) {
    for(t in seq_len(T)) {
      forward_vars[[t]] <- exp(forward_vars[[t]])
    }
  }
  
  return(forward_vars)
}

calc_HMM_forward_backward_vars_given_log_obs_probs <- function(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = TRUE) {
  # S is the number of elements in the state space
  # initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
  # trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
  
  state_space <- matrix(seq_len(S))
  
  forward_vars <- calc_HMM_forward_vars_given_log_obs_probs(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = TRUE)
  
  backward_vars <- vector("list", T)
  
  backward_vars[[T]] <- matrix(0, nrow = 1, ncol = S)
  for(t in rev(1:(T-1)) ) {
    backward_vars[[t]] <- backward_vars[[t+1]]
    for(r in 1:S) {
      temp <- backward_vars[[t+1]] + log_trans_matrix[r, ] + log_obs_probs[, t+1]
      dim(temp) <- c(1, S)
      
      backward_vars[[t]][1, r] <- logspace_sum_matrix_rows(temp)
    }
  }
  
  if(!log) {
    for(t in seq_len(T)) {
      forward_vars[[t]] <- exp(forward_vars[[t]])
      backward_vars[[t]] <- exp(backward_vars[[t]])
    }
  }
  
  return(list(forward_vars = forward_vars, backward_vars = backward_vars))
}

