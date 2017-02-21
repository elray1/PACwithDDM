rm(list = ls())
options(error = recover)

## Arguments

args <- commandArgs(trailingOnly = TRUE)

sim_pars <- list(sim_ind = as.integer(args[1]),
	obs_dist_normal = as.logical(args[2]),
	time_dependence = as.logical(args[3]),
	S = 3L,
	D = 50L,
	N_train = 50L,
	T_train = rep(200L, 50L),
	fit_method = args[4])

sim_pars <- list(sim_ind = 1L,
	obs_dist_normal = TRUE,
	time_dependence = TRUE,
	S = 3L,
	D = 50L,
	N_train = 50L,
	T_train = rep(200L, 50L),
	fit_method = "normalFMM")


## Load libraries and initialize cluster with snowfall

library("snowfall")

if(sim_pars$fit_method %in% c("parametricBoostCRF", "parametricBoostMLR")) {
	sfInit(parallel = TRUE, cpus = 16, type = "SOCK")
} else {
	sfInit(parallel = FALSE, cpus = 1, type = "SOCK")
}

sfLibrary("plyr", character.only = TRUE)

sfLibrary("MCMCpack", character.only = TRUE)
sfLibrary("mvtnorm", character.only = TRUE)
sfLibrary("tmvtnorm", character.only = TRUE)
sfLibrary("rstream", character.only = TRUE)

sfLibrary("randomForest", character.only = TRUE)
#sfLibrary("rpart", character.only = TRUE)
sfLibrary("nnet", character.only = TRUE)

if(sim_pars$fit_method %in% c("normalFMM", "normalHMM")) {
	library("mclust")
	library("car")
	library("mvtnorm")
}

sfLibrary("PACwithDDM", character.only = TRUE)

source(file.path(find.package("PACwithDDM"), "simStudyScripts", "run", "simStudyObsDist.R"))

## save path for results
if(sim_pars$obs_dist_normal) {
	obs_dist_folder <- "obsDistNormal"
} else {
	obs_dist_folder <- "obsDistNonNormal"
}

if(sim_pars$time_dependence) {
	time_dependence_folder <- "timeDependence"
} else {
	time_dependence_folder <- "noTimeDependence"
}

save_file_name <- paste0("sim", sim_pars$sim_ind, ".Rdata")


save_path <- file.path("C:", "Stat", "HMM", "PACwithDDM", "inst", "results", "SimStudy", obs_dist_folder, time_dependence_folder, sim_pars$fit_method, save_file_name)
#save_path <- file.path("/home", "er71a", "simStudy", "results", obs_dist_folder, time_dependence_folder, sim_pars$fit_method, save_file_name)
#save_path <- file.path("/home", "em"/", "Stat", "hmm", "hmmensembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, sim_pars$fit_method, save_file_name)


## Initialize rng and generate data

# load rng stream for data generation
stream_filename <- paste("rngstream_simStudy", obs_dist_folder, time_dependence_folder,
	"sim", sim_pars$sim_ind, sep = "_")

load(file = file.path(find.package("PACwithDDM"), "simStudyScripts", "rngstreams", "data_gen",
  paste0(stream_filename, ".rdata")))

# set rstream as the RNG used by R, using the rngstream object
rstream.packed(rngstream) <- FALSE
rstream.RNG(rngstream)

# initial state distribution and transition matrix for data simulation
if(sim_pars$time_dependence) {
  if(sim_pars$S == 3) {
    init_state_probs_train <- rep(1 / sim_pars$S, sim_pars$S)
    transition_matrix_train <- matrix(0.1, nrow = sim_pars$S, ncol = sim_pars$S)
    diag(transition_matrix_train) <- 0.8
  } else {
    stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
  }
} else {
  if(sim_pars$S == 3) {
    init_state_probs_train <- rep(1 / sim_pars$S, sim_pars$S)
    transition_matrix_train <- matrix(1 / sim_pars$S, nrow = sim_pars$S, ncol = sim_pars$S)
  } else {
    stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
  }
}

init_state_probs_test <- init_state_probs_train
transition_matrix_test <- transition_matrix_train

# observation distribution parameters for data simulation
obs_dist_params_train <- get_obs_dist_params(obs_dist_normal = sim_pars$obs_dist_normal,
	redundant_features_informative = TRUE, bayes_error_rate_high = FALSE)
obs_dist_params_train$D <- sim_pars$D
obs_dist_params_test <- obs_dist_params_train

# simulate training and evaluation data
train_data <- sim_data_from_HMM(sim_pars$N_train, sim_pars$T_train, init_state_probs_train, transition_matrix_train, obs_dist_params_train, robs_dist)
test_data <- sim_data_from_HMM(50L, rep(200L, 50), init_state_probs_test, transition_matrix_test, obs_dist_params_test, robs_dist)

# load rng stream for fit method
if(!identical(sim_pars$fit_method, "normalFMM")) {
  load(file = file.path(find.package("PACwithDDM"), "simStudyScripts", "rngstreams", sim_pars$fit_method,
    paste0(stream_filename, ".rdata")))
  
  # set rstream as the RNG used by R, using the rngstream object
  rstream.packed(rngstream) <- FALSE
  rstream.RNG(rngstream)
}


## Get model fit

# set up arguments to control the fit
rng_seed <- NULL
rstream_substream_offset <- 0

if(identical(sim_pars$fit_method, "parametricBoostCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-boost", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "sequence", M_bag = 1000, sequence_bag_sample_size = sim_pars$N_train, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = TRUE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
		beta_initialization = "MLR", beta_penalty_factor = 0,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(identical(sim_pars$fit_method, "parametricBoostMLR")) {
	crf_control <- lccrf_control(fit_method = "parametric-boost-MLR", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "sequence", M_bag = 1000, sequence_bag_sample_size = sim_pars$N_train, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = FALSE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
		beta_initialization = "MLR", beta_penalty_factor = 0,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(!(sim_pars$fit_method %in% c("RF", "normalHMM", "normalFMM"))) {
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
} else if(identical(sim_pars$fit_method, "normalFMM")) {
		fmm_test_data <- list(X = rbind.fill.matrix(lapply(test_data, function(comp) comp$X)),
				y = unlist(lapply(test_data, function(comp) comp$y)),
				subject = unlist(lapply(seq_along(test_data), function(compind) rep(compind, length(test_data[[compind]]$y)))))
		
		fmm_train_data <- list(X = rbind.fill.matrix(lapply(train_data, function(comp) comp$X)),
				y = unlist(lapply(train_data, function(comp) comp$y)),
				subject = unlist(lapply(seq_along(train_data), function(compind) rep(compind, length(train_data[[compind]]$y)))))
		
		classes <- levels(fmm_train_data$y)
		
		fit_time <- system.time({
					# perform box-cox transformation -- separate transformation for each covariate and class,
					# estimating transform parameters from training data and applying to both train and test data
					# first, center
					cv <- apply(fmm_train_data$X, 2, mean)
					fmm_train_data$X <- fmm_train_data$X - matrix(rep(cv, each = nrow(fmm_train_data$X)), nrow = nrow(fmm_train_data$X))
					fmm_test_data$X <- fmm_test_data$X - matrix(rep(cv, each = nrow(fmm_test_data$X)), nrow = nrow(fmm_test_data$X))
					
					transform_params <- lapply(seq_len(ncol(fmm_train_data$X)), function(colind) powerTransform(fmm_train_data$X[, colind] ~ fmm_train_data$y, family = "yjPower"))
					for(colind in seq_len(ncol(fmm_train_data$X))) {
						fmm_train_data$X[, colind] <- yjPower(fmm_train_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
						fmm_test_data$X[, colind] <- yjPower(fmm_test_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
					}
					
					# get estimated transition matrix
					log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(fmm_train_data$y,
							fmm_train_data$subject, length(classes), log = TRUE)
					
					# get estimated normal components for each class
					class_norms <- lapply(classes, function(class_name) {
								mclust_fit <- densityMclust(fmm_train_data$X[fmm_train_data$y == class_name, ], G = 1:9)
								M_mclust <- ncol(mclust_fit$parameters$mean)
								return(list(rho = mclust_fit$parameters$pro,
												mus = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$mean[, ind]),
												Sigmas = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$variance$sigma[, , ind])))
							})
					
					log_obs_probs <- as.matrix(as.data.frame(lapply(seq_along(classes), function(s) {
														dGMM(fmm_test_data$X, rhos = class_norms[[s]]$rho, mus = class_norms[[s]]$mus, Sigmas = class_norms[[s]]$Sigmas, log = TRUE)
													})))
				})
		
		log_class_probs <- log_obs_probs
		y_pred <- apply(log_obs_probs, 1, which.max)
	} else {
	crftrain_data <- list(X = rbind.fill.matrix(lapply(train_data, function(comp) comp$X)),
		y = unlist(lapply(train_data, function(comp) comp$y)),
		subject = unlist(lapply(seq_along(train_data), function(compind) rep(compind, length(train_data[[compind]]$y)))))

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

save(fit_time, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
