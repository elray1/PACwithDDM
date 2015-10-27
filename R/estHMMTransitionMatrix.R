estimate_transition_matrix_reduced_parameterization_states_observed <- function(y, subject, S, log = TRUE) {
  n_trans_total <- (length(y) - length(unique(subject)))
  n_trans_same <- 0
  for(t in seq_len(length(y) - 1))
    if(subject[t] == subject[t + 1])
      if(y[t] == y[t + 1])
        n_trans_same <- n_trans_same + 1
  
  trans_matrix <- matrix((n_trans_total - n_trans_same) / n_trans_total / (S - 1), nrow = S, ncol = S)
  diag(trans_matrix) <- n_trans_same/n_trans_total
  
  if(log) {
    return(log(trans_matrix))
  } else {
    return(trans_matrix)
  }
}

estimate_transition_matrix_diag_plus_stationary_parameterization_states_observed <- function(y, subject, S, log_stationary_class_probs, log = TRUE) {
  n_trans_total <- (length(y) - length(unique(subject)))
  n_trans_same <- 0
  for(t in seq_len(length(y) - 1))
    if(subject[t] == subject[t + 1])
      if(y[t] == y[t + 1])
        n_trans_same <- n_trans_same + 1
  
  log_trans_matrix <- matrix(rep(log_stationary_class_probs, S), byrow = TRUE, nrow = S, ncol = S)
  log_ks <- sapply(log_stationary_class_probs, function(pr) logspace_sub(0, pr))
  log_trans_matrix <- sweep(log_trans_matrix, 2, log_ks, `-`)

  log_C <- log(n_trans_same) - log(n_trans_total)
  log_trans_matrix <- logspace_sub(0, log_C) + log_trans_matrix
  
  diag(log_trans_matrix) <- log_C
  
  log_k2s <- logspace_sum_matrix_rows(log_trans_matrix)
  log_trans_matrix <- sweep(log_trans_matrix, 1, log_k2s, `-`)

#  probs <- exp(log_stationary_class_probs)
#  constr_mat <- matrix(NA, nrow = 1, ncol = S^2)
#  for(s in seq_len(S)) {
#    new_vec <- rep(0, S^2)
#    new_vec[(s - 1) * S + seq_len(S)] <- probs
#    new_vec[(s - 1) * S + s] <- 0
#    constr_mat <- rbind(constr_mat, new_vec)
#  }
#
#  for(s in seq_len(S)) {
#    new_vec <- rep(0, S^2)
#    new_vec[seq_len(S) + S * (seq_len(S) - 1)] <- 1
#    new_vec[s + S * (s - 1)] <- 0
#    constr_mat <- rbind(constr_mat, new_vec)
#  }
#
#  constr_mat <- constr_mat[-1, ]
  
  if(log) {
    return(log_trans_matrix)
  } else {
    return(exp(log_trans_matrix))
  }
}


estimate_transition_matrix_full_parameterization_states_observed <- function(y, subject, S, log = TRUE) {
  trans_counts <- matrix(0, nrow = S, ncol = S)
  for(t in seq_len(length(y) - 1))
    if(subject[t] == subject[t + 1])
      trans_counts[y[t], y[t + 1]] <- trans_counts[y[t], y[t + 1]] + 1
  
  row_sums <- apply(trans_counts, 1, sum)
  trans_matrix <- sweep(trans_counts, 1, row_sums, `/`)
  
  # if there are no observations for a given row, set estimated transition probabilities to
  # p(return to same state) = 0.9, other transitions equally likely
  trans_matrix[row_sums == 0, ] <- 0.1 / (S - 1)
  trans_matrix[row_sums == 0, row_sums == 0] <- 0.9
  
  if(log) {
    return(log(trans_matrix))
  } else {
    return(trans_matrix)
  }
}

