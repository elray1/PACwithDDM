rm(list = ls())
#options(error = recover)

## Arguments

args <- commandArgs(trailingOnly = TRUE)

sim_pars <- list(sim_ind = as.integer(args[1]),
	obs_dist_normal = as.logical(args[2]),
	redundant_features_informative = as.logical(args[3]),
	bayes_error_rate_high = as.logical(args[4]),
	S = 3L,
	D = 50L,
	N_train = 50L,
	T_train = rep(200L, 50L),
	fit_method = args[5])

#sim_pars <- list(sim_ind = as.integer(2),
#	obs_dist_normal = as.logical(FALSE),
#	redundant_features_informative = as.logical(FALSE),
#	bayes_error_rate_high = as.logical(FALSE),
#	S = 3L,
#	D = 50L,
#	N_train = 50L,
#	T_train = rep(200L, 50L),
#	fit_method = "normalHMM")

## Load libraries and initialize cluster with snowfall

library("snowfall")

if(identical(sim_pars$fit_method, "parametricBoostCRF")) {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
} else if(identical(sim_pars$fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
} else if(sim_pars$fit_method %in% c("RF", "RFHMM", "MLRHMM", "normalHMM", "unregularizedParametricCRF", "BayesRuleTestParams")) {
	sfInit(parallel = FALSE, cpus = 1, type = "SOCK")
} else {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
}

sfLibrary("plyr", character.only = TRUE)

sfLibrary("MCMCpack", character.only = TRUE)
sfLibrary("mvtnorm", character.only = TRUE)
sfLibrary("tmvtnorm", character.only = TRUE)
sfLibrary("rstream", character.only = TRUE)

sfLibrary("randomForest", character.only = TRUE)
sfLibrary("rpart", character.only = TRUE)
sfLibrary("nnet", character.only = TRUE)

if(identical(sim_pars$fit_method, "normalHMM")) {
	library("mclust")
	library("car")
}

sfLibrary("lccrf", character.only = TRUE)
sfLibrary("rayThesis", character.only = TRUE)



source(file.path(find.package("rayThesis"), "simStudyScripts", "run", "simStudyObsDist.R"))




## save path for results
if(sim_pars$obs_dist_normal) {
	obs_dist_folder <- "obsDistNormal"
} else {
	obs_dist_folder <- "obsDistNonNormal"
}

if(sim_pars$bayes_error_rate_high) {
	bayes_error_folder <- "BayesErrorLarge"
} else {
	bayes_error_folder <- "BayesErrorSmall"
}

if(sim_pars$redundant_features_informative) {
	redundant_features_folder <- "redundantFeaturesInformative"
} else {
	redundant_features_folder <- "redundantFeaturesNonInformative"
}

save_file_name <- paste0("sim", sim_pars$sim_ind, ".Rdata")


save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, sim_pars$fit_method, save_file_name)
#save_path <- file.path("/home", "er71a", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, sim_pars$fit_method, save_file_name)
#save_path <- file.path("/home", "em"/", "Stat", "hmm", "hmmensembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, sim_pars$fit_method, save_file_name)


## Initialize rng and generate data

# load rng stream for data generation
stream_filename <- paste("rngstream_simStudy", obs_dist_folder, bayes_error_folder, redundant_features_folder,
	"sim", sim_pars$sim_ind, sep = "_")

load(file = file.path(find.package("rayThesis"), "simStudyScripts", "rngstreams", "data_gen", paste0(stream_filename, ".rdata")))

# set rstream as the RNG used by R, using the rngstream object
rstream.packed(rngstream) <- FALSE
rstream.RNG(rngstream)

# initial state distribution and transition matrix for data simulation
#if(identical(sim_pars$dependence_type, "1HMM")) {
#	if(sim_pars$marginal_class_probs_equal) {
		if(sim_pars$S == 3) {
			init_state_probs_train <- rep(1 / sim_pars$S, sim_pars$S)
			transition_matrix_train <- matrix(0.1, nrow = sim_pars$S, ncol = sim_pars$S)
			diag(transition_matrix_train) <- 0.8
		} else {
			stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
		}
#	} else {
#		if(sim_pars$S == 3) {
#			init_state_probs_train <- c(2/3, 1/4, 1/12)
#			transition_matrix_train <- matrix(c(9/10, 11/60, 1/4, 11/160, 4/5, 1/20, 5/160, 1/60, 7/10), nrow = sim_pars$S, ncol = sim_pars$S)
#		} else {
#			stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
#		}
#	}
#} else if(identical(sim_pars$dependence_type, "2HMM")) {
#	if(sim_pars$marginal_class_probs_equal) {
#		if(sim_pars$S == 3) {
#			init_state_probs_train <- rep(1 / sim_pars$S, sim_pars$S)
#			transition_matrix_train <- matrix(0.1, nrow = sim_pars$S, ncol = sim_pars$S)
#			diag(transition_matrix_train) <- 0.8
#		} else {
#			stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
#		}
#	} else {
#		if(sim_pars$S == 3) {
#			init_state_probs_train <- c(2/3, 1/4, 1/12)
#			hmm1_transition_matrix_train <- matrix(c(9/10, 11/60, 1/4, 11/160, 4/5, 1/20, 5/160, 1/60, 7/10), nrow = 3, ncol = 3)
#
#			distribution_factors <- list(c(0, 1/2, 1/2), c(3/4, 1/4, 0), c(1/3, 1/3, 1/3),
#				c(9/100, 9/10, 1/100), c(1/6, 1/6, 2/3), c(1/3, 0, 2/3),
#				c(1/10, 1/2, 4/10), c(1/2, 1/4, 1/4), c(1/10, 8/10, 1/10))
#
#			distributed_transition_probs <- unlist(lapply(seq_along(hmm1_transition_matrix_train), function(ind) distribution_factors[[ind]] * hmm1_transition_matrix_train[ind]))
#
#			transition_matrix_train <- array(distributed_transition_probs, dim = c(3, 3, 3))
#		} else {
#			stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
#		}
#	}
#} else {
#	if(sim_pars$S == 3) {
#		init_state_probs_train <- rep(1 / 3, 3)
#		transition_matrix_train <- matrix(c(0, 3/4, 1/3, 1/2, 0, 2/3, 1/2, 1/4, 0), nrow = 3, ncol = 3)
#		duration_dist_params <- list(list(start_state = 1, end_state = 2, alpha = 0, beta = 1, r = 442413, M = 10, q = 1, s = 0),
#			list(start_state = 1, end_state = 3, alpha = 0, beta = 1.65, r = 0.45, M = 10, q = 1, s = 0),
#			list(start_state = 2, end_state = 1, alpha = 1808, beta = 148.41, r = 33.12, M = 10, q = 1, s = 0),
#			list(start_state = 2, end_state = 3, alpha = 0, beta = 7.39, r = 7.39, M = 10, q = 0.5, s = 0.69),
#			list(start_state = 3, end_state = 1, alpha = 0, beta = 22026, r = 0.61, M = 10, q = 0.62, s = 0.9),
#			list(start_state = 3, end_state = 2, alpha = 0, beta = 1, r = 442413, M = 10, q = 0.5, s = 0.9))
#	} else {
#		stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
#	}
#}

#if(sim_pars$marginal_class_probs_differ_train_test) {
#	init_state_probs_test <- rev(init_state_probs_train)
#	transition_matrix_test <- matrix(rev(transition_matrix_train), nrow = sim_pars$S, ncol = sim_pars$S)
#} else {
	init_state_probs_test <- init_state_probs_train
	transition_matrix_test <- transition_matrix_train
#}

# observation distribution parameters for data simulation
obs_dist_params_train <- get_obs_dist_params(obs_dist_normal = sim_pars$obs_dist_normal,
	redundant_features_informative = sim_pars$redundant_features_informative, bayes_error_rate_high = sim_pars$bayes_error_rate_high)
obs_dist_params_train$D <- sim_pars$D
obs_dist_params_test <- obs_dist_params_train

# simulate training and evaluation data
#if(identical(sim_pars$dependence_type, "1HMM")) {
	train_data <- sim_data_from_HMM(sim_pars$N_train, sim_pars$T_train, init_state_probs_train, transition_matrix_train, obs_dist_params_train, robs_dist)
	test_data <- sim_data_from_HMM(50L, rep(200L, 50), init_state_probs_test, transition_matrix_test, obs_dist_params_test, robs_dist)
#} else if(identical(sim_pars$dependence_type, "2HMM")) {
#	train_data <- sim_data_from_2HMM(sim_pars$N_train, sim_pars$T_train, init_state_probs_train, transition_matrix_train, obs_dist_params_train, robs_dist)
#	test_data <- sim_data_from_2HMM(50L, rep(200L, 50), init_state_probs_test, transition_matrix_test, obs_dist_params_test, robs_dist)
#} else {
#	train_data <- sim_data_from_TDGMM(sim_pars$N_train, sim_pars$T_train, init_state_probs_train, transition_matrix_train, duration_dist_params, obs_dist_params_train, robs_dist)
#	test_data <- sim_data_from_TDGMM(50L, rep(200L, 50), init_state_probs_test, transition_matrix_test, duration_dist_params, obs_dist_params_test, robs_dist)
#}

# load rng stream for fit method
if(! sim_pars$fit_method %in% c("BayesRuleTestParams", "BayesRuleTrainParams")) {
	load(file = file.path(find.package("rayThesis"), "simStudyScripts", "rngstreams", sim_pars$fit_method, paste0(stream_filename, ".rdata")))

	# set rstream as the RNG used by R, using the rngstream object
	rstream.packed(rngstream) <- FALSE
	rstream.RNG(rngstream)
}


## Get model fit

# set up arguments to control the fit
rng_seed <- NULL
rstream_substream_offset <- 0

if(identical(sim_pars$fit_method, "RFCRF")) {
	crf_control <- lccrf_control(fit_method = "rf", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = c("sequence", "timepoint"), M_bag = 1000, timepoint_bag_sample_size_proportion = 0.5, sequence_bag_sample_size = 2 * (sim_pars$N_train), sequence_bag_sample_replace = TRUE,
		update_log_rhos = TRUE, update_transition_matrix = TRUE, max_rf_refit_iter = 10, active_var_search_method = "random", max_attempts = 100,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(sim_pars$fit_method, "RFCRFsingletimepointbag")) {
	crf_control <- lccrf_control(fit_method = "rf", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = c("sequence", "timepoint"), M_bag = 1000, timepoint_bag_sample_size = sum(sim_pars$T_train), sequence_bag_sample_size = sum(sim_pars$T_train), sequence_bag_sample_replace = TRUE,
		update_log_rhos = TRUE, update_transition_matrix = TRUE, max_rf_refit_iter = 10, active_var_search_method = "random", max_attempts = 100,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(sim_pars$fit_method, "RFCRFseqbag")) {
	crf_control <- lccrf_control(fit_method = "rf", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = c("sequence"), M_bag = 1000, sequence_bag_sample_size = sim_pars$N_train, sequence_bag_sample_replace = TRUE,
		update_log_rhos = TRUE, update_transition_matrix = TRUE, max_rf_refit_iter = 10, active_var_search_method = "random", max_attempts = 100,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(identical(sim_pars$fit_method, "parametricBoostCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-boost", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "sequence", M_bag = 100, sequence_bag_sample_size = sim_pars$N_train, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = TRUE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
		beta_initialization = "MLR", beta_penalty_factor = 0,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(identical(sim_pars$fit_method, "L2RegularizedCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-L2-penalized-MLE", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "none",
		update_transition_matrix = TRUE, M_boost = 1, M_boost_search_threshold = 0, num_active_vars = sim_pars$D, active_var_search_method = "random", max_attempts = 1,
		beta_initialization = "MLR", beta_penalty_factor = "crossval-select", 
		K_crossval = 10,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(sim_pars$fit_method, "unregularizedParametricCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-L2-penalized-MLE", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "none",
		update_transition_matrix = TRUE, M_boost = 1, M_boost_search_threshold = 0, num_active_vars = sim_pars$D, active_var_search_method = "random", max_attempts = 1,
		beta_initialization = "MLR", beta_penalty_factor = 0,
		K_crossval = 1,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(sim_pars$fit_method, "gradientTreeBoostCRF")) {
	crf_control = lccrf_control(fit_method = "gradient-tree-boost", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "none",
		update_transition_matrix = TRUE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = sim_pars$D, active_var_search_method = "random", max_attempts = 1,
		max_tree_depth = "crossval-select", optim_method = "L-BFGS-B", crossval_method = "fit-all", K_crossval = 10, parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(sim_pars$fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
	crf_control <- lccrf_control(fit_method = "gradient-tree-boost", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "sequence", M_bag = 100, sequence_bag_sample_size = sim_pars$N_train, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = TRUE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 3, active_var_search_method = "random", max_attempts = 100,
		max_tree_depth = "crossval-select", optim_method = "L-BFGS-B", crossval_method = "one_set", K_crossval = 10, parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(!(sim_pars$fit_method %in% c("RF", "RFHMM", "MLRHMM", "normalHMM", "BayesRuleTestParams", "BayesRuleTrainParams"))) {
	stop("Invalid fit_method")
}


if(identical(sim_pars$fit_method, "RF")) {
	rftest_data <- data.frame(rbind.fill.matrix(lapply(test_data, function(comp) comp$X)))
	rftest_data$y <- unlist(lapply(test_data, function(comp) comp$y))

	rftrain_data <- data.frame(rbind.fill.matrix(lapply(train_data, function(comp) comp$X)))
	rftrain_data$y <- unlist(lapply(train_data, function(comp) comp$y))
	
	rstream.RNG(rngstream)

	fit_time <- system.time({
		fit <- randomForest(y ~ ., data = rftrain_data)
	})

	y_pred <- predict(fit, rftest_data)
	log_class_probs <- log(predict(fit, rftest_data, type = "prob"))
} else if(identical(sim_pars$fit_method, "RFHMM")) {
	rf_save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, "RF", save_file_name)
	rf_env <- new.env()
	load(rf_save_path, envir = rf_env)
	rf_log_class_probs <- rf_env$log_class_probs
	rf_log_class_probs <- lapply(seq_len(50), function(seqind) {
		rf_log_class_probs[(seqind - 1) * 200 + seq_len(200), ]
	})

	concat_X <- rbind.fill.matrix(lapply(train_data, function(comp) comp$X))
	concat_y <- unlist(lapply(train_data, function(comp) comp$y))
	concat_subject <- unlist(lapply(train_data, function(comp) rep(comp$subject, length(comp$y))))
	log_marginal_class_probs <- sapply(levels(concat_y), function(lvl) {log(sum(concat_y == lvl)) - log(sum(sim_pars$T_train))})
	log_trans_mat <- estimate_transition_matrix_full_parameterization_states_observed(concat_y, concat_subject, sim_pars$S, log = TRUE)
	
	new_HMM <- initialize_new_HMM_for_IACL(concat_X, sim_pars$T_train, sim_pars$S, root_lower_endpoints = NULL, root_upper_endpoints = NULL,
                                       endpoint_expansion_factor = 1, L_max = 1, N_min = 1,
                                       split_var_method = "best", K_var = 1, split_pt_method = "unif_data_based", K_pt = 1, split_eps = 0.0000001,
                                       return_method = "Xptr", random_rotation = FALSE, pre_sphere = FALSE)

	HMM_Xptr <- new_HMM$HMM

	set_HMM_pi_from_log(HMM_Xptr, log_marginal_class_probs)
	set_HMM_trans_matrix_from_log(HMM_Xptr, log_trans_mat)
	
	mcshane_fit <- list(log_marginal_class_probs = log_marginal_class_probs, log_trans_mat = log_trans_mat)
	
	mcshane_class_probs <- calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(HMM_Xptr, log_marginal_class_probs, 50, rep(200L, 50), rf_log_class_probs)
	log_class_probs <- mcshane_class_probs
	y_pred <- unlist(lapply(mcshane_class_probs, function(comp) apply(comp, 1, which.max)))
} else if(identical(sim_pars$fit_method, "MLRHMM")) {
	mlr_test_data <- list(X = rbind.fill.matrix(lapply(test_data, function(comp) comp$X)),
		y = unlist(lapply(test_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(test_data), function(compind) rep(compind, length(test_data[[compind]]$y)))))

	mlr_train_data <- list(X = rbind.fill.matrix(lapply(train_data, function(comp) comp$X)),
		y = unlist(lapply(train_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(train_data), function(compind) rep(compind, length(train_data[[compind]]$y)))))

	classes <- levels(mlr_train_data$y)

	fit_time <- system.time({
		mlr_fit <- fit_MLR(mlr_train_data$y, mlr_train_data$X)

		mlr_class_probs <- calc_classprobs_MLR(mlr_test_data$X, mlr_fit, log = FALSE)
		mlr_log_class_probs <- list(log(mlr_class_probs))

		log_marginal_class_probs <- sapply(levels(mlr_train_data$y), function(lvl) {log(sum(mlr_train_data$y == lvl)) - log(length(mlr_train_data$y))})
		# get estimated transition matrix
		log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(mlr_train_data$y,
			mlr_train_data$subject, length(classes), log = TRUE)
	
		new_HMM <- initialize_new_HMM_for_IACL(mlr_train_data$X, table(mlr_train_data$subject), length(classes), root_lower_endpoints = NULL, root_upper_endpoints = NULL,
										   endpoint_expansion_factor = 1, L_max = 1, N_min = 1,
										   split_var_method = "best", K_var = 1, split_pt_method = "unif_data_based", K_pt = 1, split_eps = 0.0000001,
										   return_method = "Xptr", random_rotation = FALSE, pre_sphere = FALSE)

		HMM_Xptr <- new_HMM$HMM

		set_HMM_pi_from_log(HMM_Xptr, log_marginal_class_probs)
		set_HMM_trans_matrix_from_log(HMM_Xptr, log_trans_matrix)
	
		mcshane_fit <- list(log_marginal_class_probs = log_marginal_class_probs, log_trans_mat = log_trans_matrix)
	})
	
	log_class_probs <- calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(HMM_Xptr, log_marginal_class_probs, 1, length(mlr_test_data$y), mlr_log_class_probs)
	
	y_pred <- unlist(lapply(log_class_probs, function(comp) apply(comp, 1, which.max)))
} else if(identical(sim_pars$fit_method, "normalHMM")) {
	hmm_test_data <- list(X = rbind.fill.matrix(lapply(test_data, function(comp) comp$X)),
		y = unlist(lapply(test_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(test_data), function(compind) rep(compind, length(test_data[[compind]]$y)))))

	hmm_train_data <- list(X = rbind.fill.matrix(lapply(train_data, function(comp) comp$X)),
		y = unlist(lapply(train_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(train_data), function(compind) rep(compind, length(train_data[[compind]]$y)))))

	classes <- levels(hmm_train_data$y)

	fit_time <- system.time({
		# perform box-cox transformation -- separate transformation for each covariate and class,
		# estimating transform parameters from training data and applying to both train and test data
		# first, center
		cv <- apply(hmm_train_data$X, 2, mean)
		hmm_train_data$X <- hmm_train_data$X - matrix(rep(cv, each = nrow(hmm_train_data$X)), nrow = nrow(hmm_train_data$X))
		hmm_test_data$X <- hmm_test_data$X - matrix(rep(cv, each = nrow(hmm_test_data$X)), nrow = nrow(hmm_test_data$X))

		transform_params <- lapply(seq_len(ncol(hmm_train_data$X)), function(colind) powerTransform(hmm_train_data$X[, colind] ~ hmm_train_data$y, family = "yjPower"))
		for(colind in seq_len(ncol(hmm_train_data$X))) {
			hmm_train_data$X[, colind] <- yjPower(hmm_train_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
			hmm_test_data$X[, colind] <- yjPower(hmm_test_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
		}

		# get estimated transition matrix
		log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(hmm_train_data$y,
			hmm_train_data$subject, length(classes), log = TRUE)

		# get estimated normal components for each class
		class_norms <- lapply(classes, function(class_name) {
			mclust_fit <- densityMclust(hmm_train_data$X[hmm_train_data$y == class_name, ], G = 1:9)
			M_mclust <- ncol(mclust_fit$parameters$mean)
			return(list(rho = mclust_fit$parameters$pro,
				mus = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$mean[, ind]),
				Sigmas = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$variance$sigma[, , ind])))
		})
  
		log_obs_probs <- t(as.matrix(as.data.frame(lapply(seq_along(classes), function(s) {
			dGMM(hmm_test_data$X, rhos = class_norms[[s]]$rho, mus = class_norms[[s]]$mus, Sigmas = class_norms[[s]]$Sigmas, log = TRUE)
		}))))
	})

	y_pred <- predict_HMM_given_log_obs_probs_multi_subject(length(classes), hmm_test_data$subject, log_obs_probs, log(table(hmm_train_data$y) / length(hmm_train_data$y)), log_trans_matrix)
	log_class_probs <- calc_HMM_class_probs_given_log_obs_probs_multi_subject(length(classes), hmm_test_data$subject, log_obs_probs, log(table(hmm_train_data$y) / length(hmm_train_data$y)), log_trans_matrix, log = TRUE)
} else if(identical(sim_pars$fit_method, "BayesRuleTestParams")) {
	concat_X <- rbind.fill.matrix(lapply(train_data, function(comp) comp$X))
	concat_y <- unlist(lapply(train_data, function(comp) comp$y))
	concat_subject <- unlist(lapply(train_data, function(comp) rep(comp$subject, length(comp$y))))

	fit_time <- system.time({
		new_HMM <- initialize_new_HMM_for_IACL(concat_X, sim_pars$T_train, sim_pars$S, root_lower_endpoints = NULL, root_upper_endpoints = NULL,
										   endpoint_expansion_factor = 1, L_max = 1, N_min = 1,
										   split_var_method = "best", K_var = 1, split_pt_method = "unif_data_based", K_pt = 1, split_eps = 0.0000001,
										   return_method = "Xptr", random_rotation = FALSE, pre_sphere = FALSE)

		HMM_Xptr <- new_HMM$HMM
	
		set_HMM_pi_from_log(HMM_Xptr, log(init_state_probs_test))
		set_HMM_trans_matrix_from_log(HMM_Xptr, log(transition_matrix_test))

		log_obs_probs_by_state <- lapply(test_data, function(data_comp) {
			return(sapply(seq_len(sim_pars$S), function(s) {
				return(apply(data_comp$X, 1, function(xrow) {dobs_dist(xrow, obs_dist_params_test[[s]], log = TRUE)}))
			}))
		})
	})

	temp <- predict_lccrf_given_obs_probs_marginal(HMM_Xptr, test_data, log_obs_probs_by_state, return_class_probs = TRUE)
	y_pred <- unlist(temp$predictions)
	log_class_probs <- temp$log_class_probs
} else {
#	crftest_data <- list(X = rbind.fill.matrix(lapply(test_data, function(comp) comp$X)),
#		y = unlist(lapply(test_data, function(comp) comp$y)),
#		subject = unlist(lapply(seq_along(test_data), function(compind) rep(compind, length(test_data[[compind]]$y)))))

	crftrain_data <- list(X = rbind.fill.matrix(lapply(train_data, function(comp) comp$X)),
		y = unlist(lapply(train_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(train_data), function(compind) rep(compind, length(train_data[[compind]]$y)))))

#	debug(rf_crf)
	fit_time <- system.time({
		crf_fit <- lccrf(data_concat = crftrain_data, crf_control = crf_control, rngstream)
	})

	temp <- predict_lccrf(test_data, crf_fit, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs = TRUE)
	y_pred <- unlist(temp$predictions)
	log_class_probs <- temp$log_class_probs
}


## Get measures of classifier quality and store results
confusion_matrix <- table(unlist(lapply(test_data, function(comp) comp$y)), y_pred)

num_correct <- sum(y_pred == unlist(lapply(test_data, function(comp) comp$y)))
prop_correct <- num_correct / length(unlist(lapply(test_data, function(comp) comp$y)))


if(identical(sim_pars$fit_method, "gradientTreeBoostCRF")) {
	crf_fit_one_component <- crf_fit
	crf_fit_one_component$K_crossval <- 1
	temp <- lapply(seq_along(crf_fit$component_fits), function(component_fit_ind) {
		crf_fit_one_component$component_fits <- crf_fit$component_fits[component_fit_ind]
		predict_lccrf(test_data, crf_fit_one_component, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs = TRUE)
	})
	
	y_pred_by_m <- lapply(temp, function(comp) { unlist(comp$predictions) })
	log_class_probs_by_m <- lapply(temp, function(comp) { comp$log_class_probs })
	
	confusion_matrix_by_m <- lapply(y_pred_by_m, function(y_pred) table(unlist(lapply(test_data, function(comp) comp$y)), y_pred))
	
	num_correct_by_m <- lapply(y_pred_by_m, function(y_pred) sum(y_pred == unlist(lapply(test_data, function(comp) comp$y))))
	prop_correct_by_m <- lapply(num_correct_by_m, function(num_correct) num_correct / length(unlist(lapply(test_data, function(comp) comp$y))))
	
	save(fit_time, y_pred_by_m, log_class_probs_by_m, confusion_matrix_by_m, num_correct_by_m, prop_correct_by_m, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
} else if(identical(sim_pars$fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
	oob_prop_correct <- crf_fit$oob_prop_correct
	save(oob_prop_correct, fit_time, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
} else if(identical(sim_pars$fit_method, "RF")) {
	save(fit, fit_time, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
} else if(identical(sim_pars$fit_method, "RFHMM")) {
	save(y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
} else {
	save(fit_time, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
}
