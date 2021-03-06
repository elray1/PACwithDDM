rm(list=ls())

## arguments
args <- commandArgs(trailingOnly = TRUE)

subj <- as.integer(args[1])
data_set <- args[2]
location <- args[3]
class_var <- args[4]
fit_method <- args[5]
#reduced_trans_mat_parameterization <- as.logical(args[6])

# subj <- 1L
 
# data_set <- "Mannini"
# location <- "ankle"
# class_var <- "y_category4"
# fit_method <- "normalFMM"
# fit_method <- "parametricBoostCRF"
# fit_method <- "parametricBoostMLR"
# reduced_trans_mat_parameterization <- FALSE

#subj <- 1L
#data_set <- "SasakiLab"
#location <- "wrist"
#class_var <- "y_category3"
#fit_method <- "L2RegularizedCRF"
#reduced_trans_mat_parameterization <- FALSE


## set up cluster
library("snowfall")

if(fit_method %in% c("normalHMM", "normalFMM", "RF")) {
	sfInit(parallel = FALSE, cpus = 1, type = "SOCK")
} else if(fit_method %in% "parametricBoostMLR") {
    sfInit(parallel = TRUE, cpus = 16, type = "SOCK")
} else {
	sfInit(parallel = TRUE, cpus = 50, type = "SOCK")
#    sfInit(parallel = FALSE, cpus = 1, type = "SOCK")
}



## load packages on cluster
sfLibrary("plyr", character.only = TRUE)

sfLibrary("PACwithDDM", character.only = TRUE)

if(identical(fit_method, "RF")) {
	sfLibrary("randomForest", character.only = TRUE)
}
#sfLibrary("rpart", character.only = TRUE)
sfLibrary("nnet", character.only = TRUE)

if(fit_method %in% c("normalHMM", "normalFMM")) {
	library("mclust")
	library("car")
	library("mvtnorm")
}

sfLibrary("MCMCpack", character.only = TRUE)
sfLibrary("rstream", character.only = TRUE)



## construct save path
save_file_name <- paste0("results_subject", subj, ".Rdata")

#save_path <- file.path("F:", "Evan", "PACwithDDM-linux", "pacwithddm", "inst", "results", data_set, location, class_var, fit_method, save_file_name)
#save_path <- file.path("C:", "Stat", "HMM", "PACwithDDM", "inst", "results", data_set, location, class_var, fit_method, save_file_name)
save_path <- file.path("/home", "er71a", "HMMapplication", "results", data_set, location, class_var, fit_method, save_file_name)
#save_path <- file.path("/home", "em"/", "Stat", "hmm", "hmmensembles", "HMMapplication", "results", data_set, location, class_var, fit_method, save_file_name)



## Get the data set specified by arguments
if(identical(location, "ankle")) {
	location_first_upper <- "Ankle"
	if(identical(data_set, "SasakiFreeLiving")) {
		fit_data <- SasakiFreeLivingAnkle
	} else if(identical(data_set, "SasakiLab")) {
		fit_data <- SasakiLabAnkle
	} else if(identical(data_set, "Mannini")) {
        fit_data <- ManniniAnkleCorrected
    } else {
		stop("Invalid data_set")
	}
} else if(identical(location, "hip")) {
	location_first_upper <- "Hip"
	if(identical(data_set, "SasakiFreeLiving")) {
		fit_data <- SasakiFreeLivingHip
	} else if(identical(data_set, "SasakiLab")) {
		fit_data <- SasakiLabHip
	} else {
		stop("Invalid data_set")
	}
} else if(identical(location, "wrist")) {
	location_first_upper <- "Wrist"
	if(identical(data_set, "SasakiFreeLiving")) {
		fit_data <- SasakiFreeLivingWrist
	} else if(identical(data_set, "SasakiLab")) {
		fit_data <- SasakiLabWrist
    } else if(identical(data_set, "Mannini")) {
        fit_data <- ManniniWristCorrected
	} else {
		stop("Invalid data_set")
	}
} else {
	stop("Invalid location")
}



## Put data in required format and split into training and test sets
N <- length(fit_data)
D <- ncol(fit_data[[1]]$X)

fit_data <- list(X = rbind.fill.matrix(lapply(fit_data, function(comp) comp$X)),
	y = unlist(lapply(fit_data, function(comp) comp[[class_var]])),
	subject = unlist(lapply(seq_along(fit_data), function(subj_ind) { rep(subj_ind, nrow(fit_data[[subj_ind]]$X)) })))


train_subjects <- unique(fit_data$subject)
train_subjects <- train_subjects[train_subjects != subj]
T_train <- as.vector(table(fit_data$subject[fit_data$subject != subj]))

train_data <- list(X = fit_data$X[fit_data$subject != subj, , drop = FALSE],
	y = factor(as.numeric(fit_data$y[fit_data$subject != subj]), levels = sort(unique(as.numeric(fit_data$y)))),
	subject = rep(seq_along(train_subjects), T_train))

test_data <- list(X = fit_data$X[fit_data$subject == subj, , drop = FALSE],
	y = factor(as.numeric(fit_data$y[fit_data$subject == subj]), levels = sort(unique(as.numeric(fit_data$y)))),
	subject = rep(1, sum(fit_data$subject == subj)))



## load rngstream object specific to data_set, location, class_var, fit_method, reduced_trans, update_trans
if(! fit_method %in% c("normalFMM")) {
	stream_filename <- paste("rngstream", data_set, location, class_var, fit_method,
		"subj", subj, sep = "_")
    
    load(file = file.path(find.package("PACwithDDM"), "appliedPAClassificationScripts", "rngstreams", paste0(stream_filename, ".rdata")))
}



## set up arguments to control the fit
rng_seed <- NULL
rstream_substream_offset <- 0

if(identical(fit_method, "parametricBoostCRF")) {
  crf_control <- lccrf_control(fit_method = "parametric-boost", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
    bag_method = "sequence", M_bag = 10000, sequence_bag_sample_size = N - 1, sequence_bag_sample_replace = TRUE,
    update_transition_matrix = TRUE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
    beta_initialization = "MLR", beta_penalty_factor = 0,
    optim_method = "L-BFGS-B", parallel_method = "snowfall",
    rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
} else if(identical(fit_method, "parametricBoostMLR")) {
  crf_control <- lccrf_control(fit_method = "parametric-boost-MLR", reduced_trans_mat_parameterization = FALSE, quadratic = FALSE,
    bag_method = "sequence", M_bag = 10000, sequence_bag_sample_size = N - 1, sequence_bag_sample_replace = TRUE,
    update_transition_matrix = FALSE, M_boost = 100, M_boost_search_threshold = 100, num_active_vars = 1, active_var_search_method = "random", max_attempts = 100,
    beta_initialization = "MLR", beta_penalty_factor = 0,
    optim_method = "L-BFGS-B", parallel_method = "snowfall",
    rng_method = "rstream", rng_seed = rng_seed, rstream_substream_offset = rstream_substream_offset, save_freq = Inf, save_path = "", save_filename_base = "")
} else if(!(fit_method %in% c("RF", "normalHMM", "normalFMM"))) {
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
	
} else if(identical(fit_method, "normalFMM")) {
	classes <- levels(train_data$y)

	fit_time <- system.time({
		# perform box-cox transformation -- separate transformation for each covariate and class,
		# estimating transform parameters from training data and applying to both train and test data
		transform_params <- lapply(seq_len(ncol(train_data$X)), function(colind) powerTransform(train_data$X[, colind] ~ train_data$y, family = "yjPower"))
		for(colind in seq_len(ncol(train_data$X))) {
			train_data$X[, colind] <- yjPower(train_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
			test_data$X[, colind] <- yjPower(test_data$X[, colind], coef(transform_params[[colind]], round = TRUE))
		}
    
		# get estimated normal components for each class
		class_norms <- lapply(classes, function(class_name) {
			mclust_fit <- densityMclust(train_data$X[train_data$y == class_name, ], G = 1:9)
			M_mclust <- ncol(mclust_fit$parameters$mean)
			return(list(rho = mclust_fit$parameters$pro,
				mus = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$mean[, ind]),
				Sigmas = lapply(seq_len(M_mclust), function(ind) mclust_fit$parameters$variance$sigma[, , ind])))
		})
  
		log_obs_probs <- as.matrix(as.data.frame(lapply(seq_along(classes), function(s) {
			dGMM(test_data$X, rhos = class_norms[[s]]$rho, mus = class_norms[[s]]$mus, Sigmas = class_norms[[s]]$Sigmas, log = TRUE)
		})))
	})

	log_class_probs <- log_obs_probs
	y_pred <- apply(log_obs_probs, 1, which.max)
	
} else {
	fit_time <- system.time({
		crf_fit <- lccrf(data_concat = train_data, crf_control = crf_control, rngstream = rngstream)
	})

	temp <- predict_lccrf(list(list(X = test_data$X)), crf_fit, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs = TRUE)
	y_pred <- temp$predictions[[1]]
	log_class_probs <- temp$log_class_probs
}

confusion_matrix <- table(test_data$y, y_pred)

num_correct <- sum(y_pred == test_data$y)
prop_correct <- num_correct / length(test_data$y)

## save results
save(subj, fit_time, y_pred, log_class_probs, confusion_matrix, num_correct, prop_correct, file = save_path)
