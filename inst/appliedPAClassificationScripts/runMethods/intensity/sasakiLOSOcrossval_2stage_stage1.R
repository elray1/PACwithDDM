rm(list=ls())

## Arguments

args <- commandArgs(trailingOnly = TRUE)

subj1 <- as.integer(args[1])
subj2 <- as.integer(args[2])
setting <- args[3]
location <- args[4]
class_var <- args[5]
fit_method <- args[6]
reduced_trans_mat_parameterization <- as.logical(args[7])
update_trans <- as.logical(args[8])

#subj <- 6
#setting <- "lab"
#location <- "ankle"
#class_var <- "y_category3"
#fit_method <- "baggedFeatureSubsetGradientTreeBoostCRF"
#fit_method <- "normalHMM"
#reduced_trans_mat_parameterization <- TRUE
#update_trans <- TRUE


library("snowfall")

if(identical(fit_method, "parametricBoostCRF")) {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
} else if(identical(fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
} else if(fit_method %in% c("normalHMM", "RF", "RFHMM", "MLRHMM", "SVM")) {
	sfInit(parallel = FALSE, cpus = 1, type = "SOCK")
} else {
	sfInit(parallel = TRUE, cpus = 10, type = "SOCK")
}

sfLibrary("plyr", character.only = TRUE)

sfLibrary("lccrf", character.only = TRUE)
sfLibrary("rayThesis", character.only = TRUE)

sfLibrary("randomForest", character.only = TRUE)
sfLibrary("rpart", character.only = TRUE)
sfLibrary("nnet", character.only = TRUE)

if(identical(fit_method, "normalHMM")) {
	library("mclust")
	library("car")
}

sfLibrary("MCMCpack", character.only = TRUE)
sfLibrary("rstream", character.only = TRUE)


if(reduced_trans_mat_parameterization) {
	save_file_name <- "ReducedTrans"
} else {
	save_file_name <- "FullTrans"
}

if(!update_trans)
	save_file_name <- paste0(save_file_name, "NoUpdate")

save_file_name <- paste0(save_file_name, "_subject", subj1, "_", subj2, ".Rdata")

save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "HMMapplication", "Sasaki", "intensity", "results", setting, "2stage", "stage1", location, class_var, fit_method, save_file_name)
#save_path <- file.path("/home", "er71a", "HMMapplication", "Sasaki", "intensity", "results", setting, "2stage", "stage1", location, class_var, fit_method, save_file_name)
#save_path <- file.path("/home", "em"/", "Stat", "hmm", "hmmensembles", "HMMapplication", "Sasaki", "intensity", "results", setting, "2stage", "stage1", location, class_var, fit_method, save_file_name)

## Preprocess the data
if(identical(location, "ankle")) {
	location_first_upper <- "Ankle"
	if(identical(setting, "freeliving")) {
		fit_data <- SasakiFreeLivingAnkle
	} else if(identical(setting, "lab")) {
		fit_data <- SasakiLabAnkle
	} else {
		stop("Invalid setting")
	}
} else if(identical(location, "hip")) {
	location_first_upper <- "Hip"
	if(identical(setting, "freeliving")) {
		fit_data <- SasakiFreeLivingHip
	} else if(identical(setting, "lab")) {
		fit_data <- SasakiLabHip
	} else {
		stop("Invalid setting")
	}
} else if(identical(location, "wrist")) {
	location_first_upper <- "Wrist"
	if(identical(setting, "freeliving")) {
		fit_data <- SasakiFreeLivingWrist
	} else if(identical(setting, "lab")) {
		fit_data <- SasakiLabWrist
	} else {
		stop("Invalid setting")
	}
} else {
	stop("Invalid location")
}

N <- length(fit_data)
D <- ncol(fit_data[[1]]$X)

fit_data <- list(X = rbind.fill.matrix(lapply(fit_data, function(comp) comp$X)),
	y = unlist(lapply(fit_data, function(comp) comp[[class_var]])),
	subject = unlist(lapply(seq_along(fit_data), function(subj_ind) { rep(subj_ind, nrow(fit_data[[subj_ind]]$X)) })))


all_subjects <- unique(fit_data$subject)

train_subjects <- all_subjects[!(all_subjects %in% c(subj1, subj2))]
T_train <- as.vector(table(fit_data$subject[!(fit_data$subject %in% c(subj1, subj2))]))

test_subjects <- all_subjects[(all_subjects %in% c(subj1, subj2))]
T_test <- as.vector(table(fit_data$subject[(fit_data$subject %in% c(subj1, subj2))]))

train_data <- list(X = fit_data$X[!(fit_data$subject %in% c(subj1, subj2)), , drop = FALSE],
	y = factor(as.numeric(fit_data$y[!(fit_data$subject %in% c(subj1, subj2))])),
	subject = rep(seq_along(train_subjects), T_train))

test_data <- list(X = fit_data$X[(fit_data$subject %in% c(subj1, subj2)), , drop = FALSE],
	y = factor(as.numeric(fit_data$y[(fit_data$subject %in% c(subj1, subj2))])),
	subject = rep(seq_along(test_subjects), T_test))




# load rngstream object specific to setting, location, class_var, fit_method, reduced_trans, update_trans
if(! fit_method %in% c("RFHMM", "MLRHMM", "SVM")) {
	stream_filename <- paste("rngstream_Sasaki", setting, location, fit_method,
		"case_2stage_stage1", "subj", subj1, subj2, sep = "_")

	load(file = file.path(find.package("rayThesis"), "appliedPAClassificationScripts", "Sasaki", "intensity", "rngstreams", paste0(stream_filename, ".rdata")))
}

rng_seed <- NULL
rstream_substream_offset <- 0

# set up arguments to control the fit
if(identical(fit_method, "RFCRF")) {
	crf_control <- lccrf_control(fit_method = "rf", reduced_trans_mat_parameterization = reduced_trans_mat_parameterization, quadratic = FALSE,
		bag_method = c("sequence", "timepoint"), M_bag = 1000, timepoint_bag_sample_size_proportion = 0.5, sequence_bag_sample_size = 2 * (N - 1), sequence_bag_sample_replace = TRUE,
		update_log_rhos = TRUE, update_transition_matrix = update_trans, max_rf_refit_iter = 10, active_var_search_method = "random", max_attempts = 100,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(fit_method, "RFCRFseqbag")) {
	crf_control <- lccrf_control(fit_method = "rf", reduced_trans_mat_parameterization = reduced_trans_mat_parameterization, quadratic = FALSE,
		bag_method = c("sequence"), M_bag = 1000, sequence_bag_sample_size = N - 1, sequence_bag_sample_replace = TRUE,
		update_log_rhos = TRUE, update_transition_matrix = update_trans, max_rf_refit_iter = 10, active_var_search_method = "random", max_attempts = 100,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(identical(fit_method, "parametricBoostCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-boost", reduced_trans_mat_parameterization = reduced_trans_mat_parameterization, quadratic = FALSE,
		bag_method = "sequence", M_bag = 100, sequence_bag_sample_size = N - 1, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = update_trans, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
		beta_initialization = "MLR", beta_penalty_factor = 0,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")

} else if(identical(fit_method, "L2RegularizedCRF")) {
	crf_control <- lccrf_control(fit_method = "parametric-L2-penalized-MLE", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
		bag_method = "none",
		update_transition_matrix = update_trans, M_boost = 1, M_boost_search_threshold = 0, num_active_vars = D, active_var_search_method = "random", max_attempts = 1,
		beta_initialization = "MLR", beta_penalty_factor = "crossval-select", 
		K_crossval = 10,
		optim_method = "L-BFGS-B", parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(fit_method, "gradientTreeBoostCRF")) {
	crf_control = lccrf_control(fit_method = "gradient-tree-boost", reduced_trans_mat_parameterization = reduced_trans_mat_parameterization, quadratic = FALSE,
		bag_method = "none",
		update_transition_matrix = update_trans, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = D, active_var_search_method = "random", max_attempts = 1,
		max_tree_depth = "crossval-select", optim_method = "L-BFGS-B", crossval_method = "fit-all", K_crossval = 10, parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(identical(fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
	crf_control <- lccrf_control(fit_method = "gradient-tree-boost", reduced_trans_mat_parameterization = reduced_trans_mat_parameterization, quadratic = FALSE,
		bag_method = "sequence", M_bag = 100, sequence_bag_sample_size = N - 1, sequence_bag_sample_replace = TRUE,
		update_transition_matrix = update_trans, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 3, active_var_search_method = "random", max_attempts = 100,
		max_tree_depth = "crossval-select", optim_method = "L-BFGS-B", crossval_method = "one_set", K_crossval = 10, parallel_method = "snowfall",
		rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
	
} else if(!(fit_method %in% c("RF", "RFHMM", "MLRHMM", "normalHMM"))) {
	stop("Invalid fit_method")
}


if(identical(fit_method, "RF")) {
	rftest_data <- as.data.frame(test_data$X)
	rftest_data$y <- test_data$y

	rftrain_data <- as.data.frame(train_data$X)
	rftrain_data$y <- train_data$y
	
	rstream.packed(rngstream) <- FALSE
	rstream.RNG(rngstream)

	fit_time <- system.time({
		fit <- randomForest(y ~ ., data = rftrain_data)
	})

	y_pred <- predict(fit, rftest_data)
	log_class_probs <- log(predict(fit, rftest_data, type = "prob"))

	log_class_probs <- lapply(unique(test_data$subject), function(test_subj) log_class_probs[test_data$subject == test_subj, ])
} else if(identical(fit_method, "RFHMM")) {
	rf_save_file_name <- paste0("ReducedTrans_subject", subj1, "_", subj2, ".Rdata")
	rf_save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "HMMapplication", "Sasaki", "intensity", "results", setting, "2stage", "stage1", location, class_var, "RF", rf_save_file_name)
	rf_env <- new.env()
	load(rf_save_path, envir = rf_env)
	rf_log_class_probs <- rf_env$log_class_probs

	classes <- levels(train_data$y)

	fit_time <- system.time({
		log_marginal_class_probs <- sapply(levels(train_data$y), function(lvl) {log(sum(train_data$y == lvl)) - log(length(train_data$y))})
		# get estimated transition matrix
		if(reduced_trans_mat_parameterization) {
			log_trans_matrix <- estimate_transition_matrix_reduced_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
		} else {
			log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
			log_trans_matrix[log_trans_matrix < -100] <- -100
		}
	
		new_HMM <- initialize_new_HMM_for_IACL(train_data$X, table(train_data$subject), length(classes), root_lower_endpoints = NULL, root_upper_endpoints = NULL,
										   endpoint_expansion_factor = 1, L_max = 1, N_min = 1,
										   split_var_method = "best", K_var = 1, split_pt_method = "unif_data_based", K_pt = 1, split_eps = 0.0000001,
										   return_method = "Xptr", random_rotation = FALSE, pre_sphere = FALSE)

		HMM_Xptr <- new_HMM$HMM

		set_HMM_pi_from_log(HMM_Xptr, log_marginal_class_probs)
		set_HMM_trans_matrix_from_log(HMM_Xptr, log_trans_matrix)
	
		mcshane_fit <- list(log_marginal_class_probs = log_marginal_class_probs, log_trans_mat = log_trans_matrix)
	
		log_class_probs <- lapply(rf_log_class_probs, function(rf_log_class_probs_one_subj) { calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(HMM_Xptr, log_marginal_class_probs, 1, nrow(rf_log_class_probs_one_subj), list(rf_log_class_probs_one_subj))[[1]] })
	})

	y_pred <- unlist(lapply(log_class_probs, function(comp) apply(comp, 1, which.max)))
} else if(identical(fit_method, "MLRHMM")) {
	classes <- levels(train_data$y)

	fit_time <- system.time({
		mlr_fit <- fit_MLR(train_data$y, train_data$X)

		mlr_class_probs <- calc_classprobs_MLR(test_data$X, mlr_fit, log = FALSE)
		mlr_log_class_probs <- list(log(mlr_class_probs))

		log_marginal_class_probs <- sapply(levels(train_data$y), function(lvl) {log(sum(train_data$y == lvl)) - log(length(train_data$y))})
		# get estimated transition matrix
		if(reduced_trans_mat_parameterization) {
			log_trans_matrix <- estimate_transition_matrix_reduced_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
		} else {
			log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
			log_trans_matrix[log_trans_matrix < -100] <- -100
		}
	
		new_HMM <- initialize_new_HMM_for_IACL(train_data$X, table(train_data$subject), length(classes), root_lower_endpoints = NULL, root_upper_endpoints = NULL,
										   endpoint_expansion_factor = 1, L_max = 1, N_min = 1,
										   split_var_method = "best", K_var = 1, split_pt_method = "unif_data_based", K_pt = 1, split_eps = 0.0000001,
										   return_method = "Xptr", random_rotation = FALSE, pre_sphere = FALSE)
	
		HMM_Xptr <- new_HMM$HMM

		set_HMM_pi_from_log(HMM_Xptr, log_marginal_class_probs)
		set_HMM_trans_matrix_from_log(HMM_Xptr, log_trans_matrix)
	
		mcshane_fit <- list(log_marginal_class_probs = log_marginal_class_probs, log_trans_mat = log_trans_matrix)
	
		log_class_probs <- calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(HMM_Xptr, log_marginal_class_probs, 1, length(test_data$y), mlr_log_class_probs)
	})

	y_pred <- unlist(lapply(log_class_probs, function(comp) apply(comp, 1, which.max)))
} else if(identical(fit_method, "normalHMM")) {
	classes <- levels(train_data$y)

	fit_time <- system.time({
		# perform box-cox transformation -- separate transformation for each covariate and class,
		# estimating transform parameters from training data and applying to both train and test data
		transform_params <- lapply(seq_len(ncol(train_data$X)), function(colind) powerTransform(train_data$X[, colind] ~ train_data$y, family = "yjPower"))
		for(colind in seq_len(ncol(train_data$X))) {
			train_data$X[, colind] <- yjPower(train_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
			test_data$X[, colind] <- yjPower(test_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
		}

		# get estimated transition matrix
		if(reduced_trans_mat_parameterization) {
			log_trans_matrix <- estimate_transition_matrix_reduced_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
		} else {
			log_trans_matrix <- estimate_transition_matrix_full_parameterization_states_observed(train_data$y,
				train_data$subject, length(classes), log = TRUE)
		}

		# get estimated normal components for each class
		class_norms <- lapply(classes, function(class_name) {
			mclust_fit <- densityMclust(train_data$X[train_data$y == class_name, ], G = 1:9)
			M_mclust <- ncol(mclust_fit$parameters$mean)
			return(list(rho = mclust_fit$parameters$pro,
				mus = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$mean[, ind]),
				Sigmas = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$variance$sigma[, , ind])))
		})
  
		log_obs_probs <- t(as.matrix(as.data.frame(lapply(seq_along(classes), function(s) {
			dGMM(test_data$X, rhos = class_norms[[s]]$rho, mus = class_norms[[s]]$mus, Sigmas = class_norms[[s]]$Sigmas, log = TRUE)
		}))))
	})

	y_pred <- predict_HMM_given_log_obs_probs_multi_subject(length(classes), test_data$subject, log_obs_probs, log(table(train_data$y) / length(train_data$y)), log_trans_matrix)
	log_class_probs <- calc_HMM_class_probs_given_log_obs_probs_multi_subject(length(classes), test_data$subject, log_obs_probs, log(table(train_data$y) / length(train_data$y)), log_trans_matrix, log = TRUE)
	log_class_probs <- lapply(c(subj1, subj2), function(subj) {
		log_class_probs[test_data$subject == subj, ]
	})
} else if(identical(fit_method, "SVM")) {
	svmtest_data <- as.data.frame(test_data$X)
	svmtest_data$y <- test_data$y

	svmtrain_data <- as.data.frame(train_data$X)
	svmtrain_data$y <- train_data$y

	fit_time <- system.time({
		fit <- svm(y ~ ., data = svmtrain_data, kernel = "radial", cost = 100, gamma = 0.1)
	})
	y_pred <- predict(fit, svmtest_data)
	log_class_probs <- NA
} else {
	fit_time <- system.time({
		crf_fit <- lccrf(data_concat = train_data, crf_control = crf_control, rngstream = rngstream)
	})

	test_data <- lapply(seq_len(2), function(subj) {
		list(X = test_data$X[test_data$subject == subj, ],
		       y = test_data$y[test_data$subject == subj])
	})

	temp <- predict_lccrf(test_data, crf_fit, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs = TRUE)

	log_class_probs <- temp$log_class_probs
}



save(subj1, subj2, log_class_probs, file = save_path)
