rm(list = ls())

library("PACwithDDM")
library("irr")
library("plyr")
library("mvtnorm")
library("tmvtnorm")
library("rstream")

options(warn = 2, error = recover)

source(file.path(find.package("PACwithDDM"), "simStudyScripts", "run", "simStudyObsDist.R"))

pacwithddm_source_location <- "C:/Stat/HMM/PACwithDDM"

combined_results_dir <- file.path(pacwithddm_source_location, "data")
results_dir_base <- file.path(pacwithddm_source_location, "inst/results/SimStudy")

results_df <- data.frame(obs_dist_normal = NA, time_dep = NA, fit_method = NA, sim_ind = NA,
	prop_correct = NA, F1_score_macro = NA, mse_pred = NA, mse_pred_rel_Bayes = NA)
row_num <- 1

sim_inds <- seq_len(50)
S <- 3L
D <- 50L

for(obs_dist_normal in c(FALSE, TRUE)) {
	for(time_dependence in c(FALSE, TRUE)) {
		## save path for results
		if(obs_dist_normal) {
			obs_dist_folder <- "obsDistNormal"
		} else {
			obs_dist_folder <- "obsDistNonNormal"
		}
	  
		if(time_dependence) {
			time_dependence_folder <- "timeDependence"
		} else {
			time_dependence_folder <- "noTimeDependence"
		}
    
		for(sim_ind in sim_inds) {
			# re-generate test data
		  sim_pars <- list(sim_ind = sim_ind,
      	obs_dist_normal = obs_dist_normal,
      	time_dependence = time_dependence,
      	S = 3L,
      	D = 50L,
      	N_train = 50L,
      	T_train = rep(200L, 50L),
      	fit_method = NA)
		  
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
      y_test <- unlist(lapply(test_data, function(comp) comp$y ))
      
      for(fit_method in c("normalFMM", "normalHMM", "RF", "parametricBoostCRF", "parametricBoostMLR")) {
        save_file_name <- paste0("sim", sim_ind, ".Rdata")
        save_path <- file.path(results_dir_base, obs_dist_folder, time_dependence_folder, fit_method, save_file_name)
        load(save_path)
        
        num_classes <- 3
        
        temp <- as.data.frame(confusion_matrix)
        inds <- as.matrix(temp[, 1:2])
        storage.mode(inds) <- "integer"
        vals <- temp[, 3]
        
        confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
        confusion_matrix[inds] <- temp[, 3]
        
        tp_by_class <- diag(confusion_matrix)
        tn_by_class <- sapply(seq_along(tp_by_class), function(ind) sum(tp_by_class[-ind]))
        fp_by_class <- apply(confusion_matrix, 2, sum) - tp_by_class
        fn_by_class <- apply(confusion_matrix, 1, sum) - tp_by_class
        
        # calculate macro F1
        to_keep <- (apply(confusion_matrix, 1, sum) > 0)
        precision_by_class <- tp_by_class[to_keep] / (tp_by_class[to_keep] + fp_by_class[to_keep])
        precision_by_class[is.na(precision_by_class)] <- 0
        recall_by_class <- tp_by_class[to_keep] / (tp_by_class[to_keep] + fn_by_class[to_keep])
        recall_by_class[is.na(recall_by_class)] <- 0
        
        precision_macro <- mean(precision_by_class)
        recall_macro <- mean(recall_by_class)
        F1_score_macro <- 2 * precision_macro * recall_macro / (precision_macro + recall_macro)
        
        # calculate MSE / Brier score
        est_log_class_probs <- rbind.fill.matrix(log_class_probs)
        est_class_probs <- exp(est_log_class_probs)
        
        y_test_as_matrix <- matrix(0, nrow = length(y_test), ncol = num_classes)
        y_test_as_matrix[cbind(seq_along(y_test), as.integer(y_test))] <- 1
        mse_pred <- sum((y_test_as_matrix - est_class_probs)^2) / length(y_test)
        
        results_df[row_num, ] <- list(
          obs_dist_normal = obs_dist_normal,
          time_dependence = time_dependence,
          fit_method = fit_method,
          sim_ind = sim_ind,
          prop_correct = prop_correct,
          F1_score_macro = F1_score_macro,
          mse_pred = mse_pred)
        row_num <- row_num + 1
      }
		}
	}
}

results_df$fit_method[results_df$fit_method == "normalFMM"] <- "FMM"
results_df$fit_method[results_df$fit_method == "normalHMM"] <- "HMM"
results_df$fit_method[results_df$fit_method == "RF"] <- "RF"
results_df$fit_method[results_df$fit_method == "parametricBoostCRF"] <- "CRF"
results_df$fit_method[results_df$fit_method == "parametricBoostMLR"] <- "MLR"
results_df$fit_method <- factor(results_df$fit_method, levels = c("CRF", "HMM", "MLR", "FMM", "RF"))

simStudyResults <- results_df
for(colind in seq_len(4))
	simStudyResults[, colind] <- as.factor(simStudyResults[, colind])

save(simStudyResults, file = file.path(combined_results_dir, "simStudyResults.rdata"))
