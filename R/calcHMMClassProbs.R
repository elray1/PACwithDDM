calc_HMM_class_probs_given_log_obs_probs <- function(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = FALSE) {
  # S is the number of elements in the state space
  # initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
  # trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
  fb <- calc_HMM_forward_backward_vars_given_log_obs_probs(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = TRUE)
  
  ps <- matrix(NA, nrow = T, ncol = S)
  
  norm_const <- logspace_sum_matrix_rows(matrix(fb$forward_vars[[T]], nrow = 1))
  
  for(t in 1:T) {
    temp <- fb$forward_vars[[t]] + fb$backward_vars[[t]]
    
    ps[t, ] <- temp - norm_const
  }
  
  if(log) {
    return(ps)
  } else {
    return(exp(ps))
  }
}

calc_HMM_class_probs_given_log_obs_probs_multi_subject <- function(S, subject_inds, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = FALSE) {
  # S is the number of elements in the state space
  # initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
  # trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
  ps <- matrix(NA, nrow = length(subject_inds), ncol = S)
  
  for(subject in unique(subject_inds)) {
    inds <- which(subject_inds == subject)
    ps[inds, ] <- calc_HMM_class_probs_given_log_obs_probs(length(inds), S, log_obs_probs[ , inds, drop = FALSE], log_initial_state_distn, log_trans_matrix, log)
  }
  
  return(ps)
}

calc_HMM_class_probs <- function(X, S, log_initial_state_distn, log_trans_matrix, dEstDist, est_params, log = FALSE) {
	# X is a T by m matrix, where row t contains the observed values of the covariates x_{t, 1}, ..., x_{t, m} for one subject
	# S is the number of elements in the state space
	# initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
	# trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
	
	T <- nrow(X)
	
	state_space <- matrix(seq_len(S))
	
	# An S by T matrix, where entry (s, t) is the likelihood of the observed data at time t if the state at the time was state s.
	log_obs_probs <- t(as.matrix(as.data.frame(lapply(state_space, function(s) {dEstDist(X, est_params[[s]], log = TRUE)}))))
#	log_obs_probs <- apply(matrix(seq_len(T)), 1, function(t) {apply(state_space, 1, function(s) {
#							dEstDist(x = X[t, , drop = FALSE], est_params[[s]], log=TRUE)
#						})})

  return(calc_HMM_class_probs_given_log_obs_probs(T, S, log_obs_probs, log_initial_state_distn, log_trans_matrix, log))
}

calc_rfdeHMM_class_probs <- function(X, S, initial_state_distn, trans_matrix, class_rfs, class_max_depths, log = FALSE) {
  # X is a T by m matrix, where row t contains the observed values of the covariates x_{t, 1}, ..., x_{t, m} for one subject
  # S is the number of elements in the state space
  # initial_state_distn is a vector of length S containing the probability distribution for the initial state of the hidden Markov chain
  # trans_matrix is an S by S matrix containing the transition probabilities for the hidden Markov chain
  
  T <- nrow(X)
  
  state_space <- matrix(seq_len(S))
  
  # An S by T matrix, where entry (s, t) is the likelihood of the observed data at time t if the state at the time was state s.
  log_obs_probs <- t(as.matrix(as.data.frame(lapply(state_space, function(s) {drfde(X, class_rfs[[s]], class_max_depths[[s]], log = TRUE)}))))
  
  return(calc_HMM_class_probs_given_log_obs_probs(T, S, log_obs_probs, initial_state_distn, trans_matrix, log))
}

predict_tedeHMM <- function(X, S, initial_state_distn, trans_matrix, class_rfs, class_max_depths) {
  class_probs <- calc_rfdeHMM_class_probs(X, S, initial_state_distn, trans_matrix, class_rfs, class_max_depths, log = TRUE)
  
  return(apply(class_probs, 1, which.max))
}

predict_tedeHMM_multi_subject <- function(X, subject, S, initial_state_distn, trans_matrix, class_rfs, class_max_depths) {
  unique_subjects <- unique(subject)
  y_pred_by_subject <- lapply(unique_subjects, function(sub) {
    sub_inds <- (subject == sub)
    predict_tedeHMM(X[sub_inds, , drop = FALSE], S, initial_state_distn, trans_matrix, class_rfs, class_max_depths)
  })
  
  return(unlist(y_pred_by_subject))
}

predict_HMM_given_log_obs_probs_multi_subject <- function(S, subject_inds, log_obs_probs, log_initial_state_distn, log_trans_matrix) {
  class_probs <- calc_HMM_class_probs_given_log_obs_probs_multi_subject(S, subject_inds, log_obs_probs, log_initial_state_distn, log_trans_matrix, log = FALSE)
  
  return(apply(class_probs, 1, which.max))
}
