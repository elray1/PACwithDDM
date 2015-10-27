calc_log_McShane_class_probs_given_log_static_class_probs_R_interface <- function(HMM_Xptr, log_marginal_class_probs, N, T_sub, log_static_class_probs) {
  return(.Call("calc_log_McShane_class_probs_given_log_static_class_probs_C",
               HMM_Xptr, as.numeric(log_marginal_class_probs), as.integer(N), as.integer(T_sub), log_static_class_probs))
}
