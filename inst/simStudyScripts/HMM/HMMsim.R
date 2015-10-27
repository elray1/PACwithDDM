options(error = recover)
rm(list = ls())

#save_location <- "~/Stat/hmm/hmmensembles/simResults/HMM/sim/"
#save_location <- "C:/Stat/HMM/HMMEnsembles/simResults/HMM/sim/"
save_location <- "/home/er71a/HMMSim/sim/results/"


library("snowfall")

args <- commandArgs(trailingOnly = TRUE)

sim_start <- as.integer(args[1])
sim_end <- as.integer(args[2])
T_train <- as.integer(args[3])
D <- as.integer(args[4])
idim <- args[5]
S <- as.integer(args[6])

sim_start <- 1L
sim_end <- 1L
T_train <- 2000L
D <- 1L
idim <- "F"
S <- 2L

n_sims <- sim_end - sim_start + 1

if(identical(D, 1L)) {
  n_cores <- 1L
} else {
  n_cores <- 2L
}
#n_cores <- 1L

sfInit(parallel = TRUE, cpus = n_cores, type = "SOCK")
#sfInit(parallel = FALSE)

# load necessary libraries on the cluster nodes
sfLibrary("rstream", character.only = TRUE)

sfLibrary("rayThesisSimStudy", character.only = TRUE)
sfLibrary("TreeHMMs", character.only = TRUE)

#debug(calc_validation_criteria_one_Lmax)

sfLibrary("MCMCpack", character.only = TRUE)
sfLibrary("mvtnorm", character.only = TRUE)
sfLibrary("tmvtnorm", character.only = TRUE)

sfLibrary("robust", character.only = TRUE)

sfLibrary("mclust", character.only = TRUE)
sfLibrary("randomForest", character.only = TRUE)



# all combinations simulation parameters
all_sim_pars <- data.frame(sim_ind = rep(seq(from = sim_start, to = sim_end, by = 1), times = 2),
                           D = rep(D, times = 2 * n_sims),
                           T_train = rep(T_train, times = 2 * n_sims),
                           S = rep(S, times = 2 * n_sims),
                           transition_matrix_parameterization = rep("full", times = 2 * n_sims),
                           idim = rep(c("F", "S"), each = n_sims),
                           stringsAsFactors = FALSE)

inds_rem <- which(all_sim_pars$D < 6 & all_sim_pars$idim == "S")
if(length(inds_rem) > 0) {
  all_sim_pars <- all_sim_pars[-inds_rem, ]
}

all_sim_pars <- lapply(seq_len(nrow(all_sim_pars)), function(ind) {
  all_sim_pars[ind, ]
})



# create a new rstream object. The seed is a vector of 6 integers
set.seed(3510)
rngstream <-new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

# pack the rstream object in preparation for exporting it to cluster nodes
rstream.packed(rngstream) <- TRUE

#results <- sfSapply(all_sim_pars, function(sim_pars) {
#results <- sapply(all_sim_pars, function(sim_pars) {
sim_pars <- all_sim_pars[[1]]
  # each sim_ind has 55 * (1 + 100) = 5555 reserved RNG substreams:
  #   - 5 possible values of S:
  #     -- S = 2
  #     -- S = 5
  #     -- three extra slots reserved for future use with additional values of S
  #   - For each value of S, 11 combinations of D and idim:
  #     -- D = 1, idim = F
  #     -- D = 6, idim = S
  #     -- D = 6, idim = F
  #     -- D = 10, idim = S
  #     -- D = 10, idim = F
  #     -- D = 30, idim = S
  #     -- D = 30, idim = F
  #     -- D = 50, idim = S
  #     -- D = 50, idim = F
  #     -- two extra slots reserved for future use with additional values of D
  #   - For each scenario (combination of S, D and idim), 1 substream is reserved for generating parameter values.
  #     note that for a given scenario, we use the same parameter values for all N_train and bdd
  #   - For each scenario (combination of S, D and idim), 100 substreams are reserved in 10 blocks of 10 substreams:
  #     -- Each block is for a unique value of T_train:
  #       --- T_train = S * 10^2; T_train = S * 10^3; T_train = S * 10^4
  #       --- 7 extra blocks in case I decide to add something else later
  #     -- Each block contains 10 substreams, for use to generate data or perform estimation with a unique density estimation method
  #       --- 1 for data generation
  #       --- 2 for RFHMM
  #       --- 3 for JboostHMM
  #       --- 4 for CboostHMM
  #       --- 5 for mclust+HMM
  #       --- 6 for RF
  #       --- 7 for McShane
  #       --- 8 - 10 reserved for future use
  
  # N_eval = size of importance samples used in estimating evaluation criteria for density estimates
  N_train <- 10
  N_eval <- 20
  T_eval <- 5000
  
  # in the cluster node where this code block is being executed,
  # unpack the exported rstream object
  rstream.packed(rngstream) <- FALSE
  
  # advance to the start of the 5555 consecutive substreams reserved for this sim_ind
  for(i in seq_len(1 + (sim_pars$sim_ind - 1) * 5555))
    rstream.nextsubstream(rngstream)
  
  # advance to the start of the set of 101 substreams reserved for this scenario
  #     -- D = 1, idim = F
  #     -- D = 6, idim = S
  #     -- D = 6, idim = F
  #     -- D = 10, idim = S
  #     -- D = 10, idim = F
  #     -- D = 30, idim = S
  #     -- D = 30, idim = F
  #     -- D = 50, idim = S
  #     -- D = 50, idim = F
  #     -- two extra slots reserved for future use with additional values of D
  # There are 11 scenarios; scenario_num is in the range of 0 to 10
  if(identical(sim_pars$S, 2L)) {
    scenario_num <- 0
  } else if(identical(sim_pars$s, 5L)) {
    scenario_num <- 11
  } else {
    stop("Invalid value for S")
  }
  
  if(identical(sim_pars$D, 1L)) {
    scenario_num <- scenario_num + 0
  } else if(identical(sim_pars$D, 6L)) {
    scenario_num <- scenario_num + 1
  } else if(identical(sim_pars$D, 10L)) {
    scenario_num <- scenario_num + 3
  } else if(identical(sim_pars$D, 30L)) {
    scenario_num <- scenario_num + 5
  } else if(identical(sim_pars$D, 50L)) {
    scenario_num <- scenario_num + 7
  } else {
    stop("Invalid value for D")
  }
  
  if(identical(sim_pars$idim, "S")) {
    scenario_num <- scenario_num + 1
  } else if(! identical(sim_pars$idim, "F")) { # note that if D == 1, sim_pars$idim always == "F"
    stop("Invalid value for idim")
  }
  
  for(i in seq_len(scenario_num * 101))
    rstream.nextsubstream(rngstream)
  
  
  # set rstream as the RNG used by R, using the rngstream object
  rstream.RNG(rngstream)
  
  # generate parameter values
  log_init_state_probs <- log(rep(1/sim_pars$S, sim_pars$S))
  
  if(identical(sim_pars$transition_matrix_parameterization, "full")) {
    log_transition_matrix <- generate_transition_matrix_full_parameterization(S = sim_pars$S, log = TRUE)
  } else if(identical(sim_pars$transition_matrix_parameterization, "reduced")) {
    log_transition_matrix <- generate_transition_matrix_reduced_parameterization(S = sim_pars$S, log = TRUE)
  } else {
    stop("Invalid value for log_transition_matrix: must be either \"full\" or \"reduced\"")
  }
  
  obs_dist_params <- list(D = sim_pars$D,
                          idim = sim_pars$idim,
                          params = r_overlapping_GMM_params(S = sim_pars$S, D = sim_pars$D, idim = sim_pars$idim, M = 5, alpha = rep(1, 5), mu_scale = 5,
                                              Xi = rep(0, sim_pars$D), Psi = diag(sim_pars$D),
                                              Sigma_scale = 0.25, nu = sim_pars$D, iwish_scale = diag(sim_pars$D),
                                              rho_adj_var = 0.1 / sim_pars$D, shift_cov_scale = sim_pars$D,
                                              cov_rotate_bounds = list(lower = -pi, upper = pi), cov_rotate_var = pi, cov_scale_var = 4,
                                              cov_expon_bounds = list(lower = 0.1, upper = 1.1), cov_expon_mean = 1, cov_expon_var = 1))

#  lapply(seq_len(2), function(s) {
#    dfn <- function(x) {dGMM(x, rho = obs_dist_params$params[[s]]$rho, mus = obs_dist_params$params[[s]]$mus, Sigmas = obs_dist_params$params[[s]]$Sigmas)}
#    curve(dfn, col = s, add = (s > 1), from = -15, to = 15)
#  })
  
  # Advance to the correct block of substreams corresponding to T_train
  #       --- T_train = S * 10^2; T_train = S * 10^3; T_train = S * 10^4
  # There are 10 blocks; block_num is in the range of 0 to 9
  if(identical(sim_pars$T_train, as.integer(sim_pars$S * 10^2))) {
    block_num <- 0
  } else if(identical(sim_pars$T_train, as.integer(sim_pars$S * 10^3))) {
    block_num <- 1
  } else if(identical(sim_pars$T_train, as.integer(sim_pars$S * 10^4))) {
    block_num <- 2
  } else {
    stop("Invalid value for T_train")
  }
  
  # advance to start of this block of substreams: 10 substreams per block
  for(i in seq_len(block_num * 10))
    rstream.nextsubstream(rngstream)
  
  # object to store results
  results <- list(sim_ind = sim_pars$sim_ind, D = sim_pars$D, T_train = sim_pars$T_train, idim = sim_pars$idim, S = sim_pars$S)
  
  
  # DATA GENERATION
  # advance to the correct substream for data generation
  rstream.nextsubstream(rngstream)
  
  # set rstream as the RNG used by R, using the rngstream object
  rstream.RNG(rngstream)
  
  # simulate training and evaluation data
  temp <- sim_data_from_HMM(N_train, sim_pars$T_train, log_init_state_probs, log_transition_matrix, obs_dist_params)
  X_train <- temp$X
  y_train <- temp$y
  subject_train <- temp$subject
  
  inds_by_state_train <- lapply(seq_len(sim_pars$S), function(s) {
    which(y_train == s)
  })
  
  temp <- sim_data_from_HMM(N_eval, T_eval, log_init_state_probs, log_transition_matrix, obs_dist_params)
  X_eval <- temp$X
  y_eval <- temp$y
  subject_eval <- temp$subject
  
  # Get a Maximum Likelihood Estimate of the transition matrix, used by all HMM-based methods below
  log_est_initial_state_probs <- rep(- log(sim_pars$S), sim_pars$S)
  log_est_transition_matrix <- estimate_transition_matrix_full_parameterization_states_observed(y_train, subject_train, sim_pars$S, log = TRUE)

  temp <- condHMMboostrfdeWithGaussianComponent(X_train, y_train, subject_train, sim_pars$S, M_overall = NULL, end_threshold = 20, min_keep_overall = 1,
                                              n_holdout = NULL, perc_holdout = NULL, K_holdout = 10, holdout_method = "averageAll",
                                              log_est_initial_state_probs, log_est_transition_matrix,
                                              M_rf = 10L, L_max = 1L, N_min = 0L,
                                              split_var_method = "best", K_var = ceiling(sim_pars$D / 2), split_pt_method = "unif_data_based", K_pt = 2,
                                              random_rotation = TRUE, bag = TRUE, pre_sphere = TRUE, shrinkage = 0.05, return_method = "Xptr", verbose = 0, lapply_fn = lapply)
  
  # RFHMM
  # advance to the correct substream for RFHMM
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
  results$train_time_RFHMM <- system.time({
    est_obs_dist_params <- vector("list", sim_pars$S)
    for(s in seq_len(sim_pars$S)) {
      est_obs_dist_params[[s]] <- crossval_select_RF_L(X_train[inds_by_state_train[[s]], , drop = FALSE], K_crossval = 10, M_xval = 200, M_final = 500, L_max_range = NULL, expand_L_max_range = TRUE,
                                                 N_min = 0, split_var_method = "best", K_var = ceiling(sim_pars$D / 2), split_pt_method = "unif_data_based", K_pt = 2,
                                                 xval_return_method = "Xptr", final_return_method = "Xptr", random_rotation = TRUE, bag = TRUE, pre_sphere = TRUE)$RF
    }
  })
  
  # Classify
  results$pred_time_RFHMM <- system.time({
    y_pred <- predict_tedeHMM_multi_subject(X_eval, subject_eval, sim_pars$S, log_est_initial_state_probs, log_est_transition_matrix, est_obs_dist_params, class_max_depths = NULL)
  })
  
  # Calculate summary statistics
  results$confusion_matrix_RFHMM <- table(y_eval, y_pred)
  
  
  
  # JboostHMM
  # advance to the correct substream for JboostHMM
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
  results$train_time_JboostHMM <- system.time({
    est_obs_dist_params <- vector("list", sim_pars$S)
    for(s in seq_len(sim_pars$S)) {
      est_obs_dist_params[[s]] <- boost_rfde(X_train[inds_by_state_train[[s]], , drop = FALSE], M_boost = 20, end_threshold = 0, min_keep = 1, n_holdout = NULL, perc_holdout = NULL, K_holdout = 5,
                                       holdout_method = "averageAll", shrinkage = 0.95, obs_weight_method = "outlier", M_rf = 10L, L_max = 1L, N_min = 0L,
                                       split_var_method = "best", K_var = ceiling(sim_pars$D / 2), split_pt_method = "unif_data_based", K_pt = 2,
                                       random_rotation = TRUE, bag = TRUE, pre_sphere = TRUE, return_method = "Xptr", verbose = 0)
    }
  })
  
  # Classify
  results$pred_time_JboostHMM <- system.time({
    y_pred <- predict_tedeHMM_multi_subject(X_eval, subject_eval, sim_pars$S, log_est_initial_state_probs, log_est_transition_matrix, est_obs_dist_params, class_max_depths = NULL)
  })
  
  # Calculate and save confusion matrix
  results$confusion_matrix_JboostHMM <- table(y_eval, y_pred)

  
  
  # CboostHMM
  # advance to the correct substream for CboostHMM
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
debug(grow_tree_C_interface)
  temp <- condHMMboostrfdeWithGaussianComponent(X_train, y_train, subject_train, sim_pars$S, M_overall = NULL, end_threshold = 20, min_keep_overall = 1,
                                              n_holdout = NULL, perc_holdout = NULL, K_holdout = 10, holdout_method = "averageAll",
                                              log_est_initial_state_probs, log_est_transition_matrix,
                                              M_rf = 10L, L_max = 1L, N_min = 0L,
                                              split_var_method = "best", K_var = ceiling(sim_pars$D / 2), split_pt_method = "unif_data_based", K_pt = 2,
                                              random_rotation = TRUE, bag = TRUE, pre_sphere = TRUE, shrinkage = 0.05, return_method = "Xptr", verbose = 0, lapply_fn = lapply)
  
  # Classify
  
  # Calculate summary statistics
  
  
  # mclust+HMM
  # advance to the correct substream for mclust + HMM
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
  results$train_time_mclustHMM <- system.time({
    mclustHMM_fit <- "blah"
  })
  
  # Classify
  results$pred_time_mclustHMM <- system.time({
    y_pred <- "blah"
  })
  
  # Calculate summary statistics
  results$confusion_matrix_mclustHMM <- table(y_eval, y_pred)
  
  # RF
  # advance to the correct substream for RF
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
  train_data_for_RF <- as.data.frame(X_train)
  names(train_data_for_RF) <- paste0("X_", seq_len(sim_pars$D))
  train_data_for_RF$y <- as.factor(y_train)
  
  results$train_time_RF <- system.time({
    rf_fit <- randomForest(y ~ ., data = train_data_for_RF)
  })
  
  # Classify
  eval_data_for_RF <- as.data.frame(X_eval)
  names(eval_data_for_RF) <- paste0("X_", seq_len(sim_pars$D))
  
  results$pred_time_RF <- system.time({
    y_pred <- predict(rf_fit, eval_data_for_RF)
  })
  
  # Calculate summary statistics
  results$confusion_matrix_RF <- table(y_eval, y_pred)
  
  
  
  # McShane
  # advance to the correct substream for McShane
  rstream.nextsubstream(rngstream)
  rstream.RNG(rngstream)
  
  # Train
  results$train_time_McShane <- system.time({
    randomForest(Y~., data=train.data)
  })
  
  # Classify
  results$pred_time_McShane <- system.time({
    y_pred <- predict(subj.fit.RF[[subj]], test.data[, !(colnames(test.data) %in% c("y", "Y"))])
  })
  
  # Calculate summary statistics
  results$confusion_matrix_McShane <- table(y_eval, y_pred)
  
  
  return(results)
#})

# stop the cluster
sfStop()
