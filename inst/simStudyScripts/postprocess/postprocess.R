rm(list = ls())

library("rayThesis")
library("irr")
library("plyr")
library("mvtnorm")
library("tmvtnorm")
library("rstream")

options(warn = 2, error = recover)

source(file.path(find.package("rayThesis"), "simStudyScripts", "run", "simStudyObsDist.R"))

results_dir <- "C:/Stat/HMM/HMMEnsembles/simStudy/results"
combined_results_dir <- "C:/Stat/HMM/HMMEnsembles/rayThesis/data/"

results_df <- data.frame(obs_dist_normal = NA, bayes_error_rate_high = NA, redundant_features_informative = NA, fit_method = NA, sim_ind = NA,
	prop_correct = NA, F1_score_macro = NA, mse_pred = NA, mse_pred_rel_Bayes = NA)
row_num <- 1

sim_inds <- seq_len(50)
S <- 3L
D <- 50L

for(obs_dist_normal in c(FALSE, TRUE)) {
	for(bayes_error_rate_high in c(FALSE, TRUE)) {
		for(redundant_features_informative in TRUE) {
#		for(redundant_features_informative in c(FALSE, TRUE)) {
			## save path for results
			if(obs_dist_normal) {
				obs_dist_folder <- "obsDistNormal"
			} else {
				obs_dist_folder <- "obsDistNonNormal"
			}

			if(bayes_error_rate_high) {
				bayes_error_folder <- "BayesErrorLarge"
			} else {
				bayes_error_folder <- "BayesErrorSmall"
			}

			if(redundant_features_informative) {
				redundant_features_folder <- "redundantFeaturesInformative"
			} else {
				redundant_features_folder <- "redundantFeaturesNonInformative"
			}

			for(sim_ind in sim_inds) {
				# re-generate test data
				# load rng stream for data generation
				stream_filename <- paste("rngstream_simStudy", obs_dist_folder, bayes_error_folder, redundant_features_folder,
					"sim", sim_ind, sep = "_")

				load(file = file.path(find.package("rayThesis"), "simStudyScripts", "rngstreams", "data_gen", paste0(stream_filename, ".rdata")))

				# set rstream as the RNG used by R, using the rngstream object
				rstream.packed(rngstream) <- FALSE
				rstream.RNG(rngstream)

				if(S == 3) {
					init_state_probs_train <- rep(1 / S, S)
					transition_matrix_train <- matrix(0.1, nrow = S, ncol = S)
					diag(transition_matrix_train) <- 0.8
				} else {
					stop("Invalid value for S: must be 3.")  # allow for the possibility of more values for S later
				}

				init_state_probs_test <- init_state_probs_train
				transition_matrix_test <- transition_matrix_train

				# observation distribution parameters for data simulation
				obs_dist_params_train <- get_obs_dist_params(obs_dist_normal = obs_dist_normal,
					redundant_features_informative = redundant_features_informative, bayes_error_rate_high = bayes_error_rate_high)
				obs_dist_params_train$D <- D
				obs_dist_params_test <- obs_dist_params_train

				# simulate training and evaluation data
				#if(identical(sim_pars$dependence_type, "1HMM")) {
				train_data <- sim_data_from_HMM(50L, rep(200L, 50), init_state_probs_train, transition_matrix_train, obs_dist_params_train, robs_dist)
				test_data <- sim_data_from_HMM(50L, rep(200L, 50), init_state_probs_test, transition_matrix_test, obs_dist_params_test, robs_dist)
				y_test <- unlist(lapply(test_data, function(comp) comp$y ))

				# load Bayes rule class probabilities
				save_file_name <- paste0("sim", sim_ind, ".Rdata")
				save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, "BayesRuleTestParams", save_file_name)
				load(save_path)
				bayes_rule_log_class_probs <- rbind.fill.matrix(log_class_probs)


				for(fit_method in c("baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "RFCRF", "RFCRFseqbag", "gradientTreeBoostCRF", "unregularizedParametricCRF", "RF", "RFHMM", "MLRHMM", "normalHMM", "BayesRuleTestParams")) {
					save_file_name <- paste0("sim", sim_ind, ".Rdata")
					save_path <- file.path("C:", "Stat", "HMM", "HMMEnsembles", "simStudy", "results", obs_dist_folder, bayes_error_folder, redundant_features_folder, fit_method, save_file_name)
					load(save_path)

					num_classes <- 3

					if(identical(fit_method, "gradientTreeBoostCRF")) {
						F1_score_macro_by_m <- lapply(confusion_matrix_by_m, function(confusion_matrix, num_classes) {
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

							return(F1_score_macro)
						}, num_classes)

						mse_pred_by_m <- lapply(log_class_probs_by_m, function(log_class_probs) {
							if(fit_method %in% c("baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "RFCRF", "RFCRFseqbag", "gradientTreeBoostCRF", "unregularizedParametricCRF", "SVM", "RFHMM", "MLRHMM")) {
								est_class_probs <- exp(rbind.fill.matrix(log_class_probs))
							} else if(fit_method %in% c("SVM")) {
								est_class_probs <- matrix(NA, nrow = length(y_pred), ncol = num_classes)
							} else {
								est_class_probs <- exp(log_class_probs)
							}
							y_test_as_matrix <- matrix(0, nrow = length(y_test), ncol = num_classes)
							y_test_as_matrix[cbind(seq_along(y_test), as.integer(y_test))] <- 1
							mse_pred <- sum((y_test_as_matrix - est_class_probs)^2) / length(y_test)

							mse_pred_rel_Bayes <- sum((est_class_probs - exp(bayes_rule_log_class_probs))^2) / length(y_test)

							return(list(mse_pred = mse_pred, mse_pred_rel_Bayes = mse_pred_rel_Bayes))
						})

						mse_pred_rel_Bayes_by_m <- lapply(mse_pred_by_m, function(comp) comp$mse_pred_rel_Bayes)
						mse_pred_by_m <- lapply(mse_pred_by_m, function(comp) comp$mse_pred)

						# worst values
						results_df[row_num, ] <- list(obs_dist_normal = obs_dist_normal, bayes_error_rate_high = bayes_error_rate_high, redundant_features_informative = redundant_features_informative, fit_method = "gradientTreeBoostCRF-worst", sim_ind = sim_ind,
							prop_correct = min(unlist(prop_correct_by_m)), F1_score_macro = min(unlist(F1_score_macro_by_m)), mse_pred = max(unlist(mse_pred_by_m)), mse_pred_rel_Bayes = max(unlist(mse_pred_rel_Bayes_by_m)))
						row_num <- row_num + 1

						# median values
						results_df[row_num, ] <- list(obs_dist_normal = obs_dist_normal, bayes_error_rate_high = bayes_error_rate_high, redundant_features_informative = redundant_features_informative, fit_method = "gradientTreeBoostCRF-med", sim_ind = sim_ind,
							prop_correct = median(unlist(prop_correct_by_m)), F1_score_macro = median(unlist(F1_score_macro_by_m)), mse_pred = median(unlist(mse_pred_by_m)), mse_pred_rel_Bayes = median(unlist(mse_pred_rel_Bayes_by_m)))
						row_num <- row_num + 1

						# best values
						results_df[row_num, ] <- list(obs_dist_normal = obs_dist_normal, bayes_error_rate_high = bayes_error_rate_high, redundant_features_informative = redundant_features_informative, fit_method = "gradientTreeBoostCRF-best", sim_ind = sim_ind,
							prop_correct = max(unlist(prop_correct_by_m)), F1_score_macro = max(unlist(F1_score_macro_by_m)), mse_pred = min(unlist(mse_pred_by_m)), mse_pred_rel_Bayes = min(unlist(mse_pred_rel_Bayes_by_m)))
						row_num <- row_num + 1
					}

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
					if(fit_method %in% c("baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "RFCRF", "RFCRFseqbag", "gradientTreeBoostCRF", "unregularizedParametricCRF", "SVM", "RFHMM", "MLRHMM")) {
						est_log_class_probs <- rbind.fill.matrix(log_class_probs)
						est_class_probs <- exp(est_log_class_probs)
					} else {
						est_log_class_probs <- rbind.fill.matrix(log_class_probs)
						est_class_probs <- exp(est_log_class_probs)
					}
					y_test_as_matrix <- matrix(0, nrow = length(y_test), ncol = num_classes)
					y_test_as_matrix[cbind(seq_along(y_test), as.integer(y_test))] <- 1
					mse_pred <- sum((y_test_as_matrix - est_class_probs)^2) / length(y_test)

					mse_pred_rel_Bayes <- sum((est_class_probs - exp(bayes_rule_log_class_probs))^2) / length(y_test)

					results_df[row_num, ] <- list(obs_dist_normal = obs_dist_normal, bayes_error_rate_high = bayes_error_rate_high, redundant_features_informative = redundant_features_informative, fit_method = fit_method, sim_ind = sim_ind,
						prop_correct = prop_correct, F1_score_macro = F1_score_macro, mse_pred = mse_pred, mse_pred_rel_Bayes = mse_pred_rel_Bayes)
					row_num <- row_num + 1
				}
			}
		}
	}
}

results_df$fit_method[results_df$fit_method == "parametricBoostCRF"] <- "BB-Par-CRF"
results_df$fit_method[results_df$fit_method == "baggedFeatureSubsetGradientTreeBoostCRF"] <- "BB-Nonpar-CRF"
results_df$fit_method[results_df$fit_method == "RFCRF"] <- "RF-CRF"
results_df$fit_method[results_df$fit_method == "RFCRFseqbag"] <- "RF-seq-CRF"
results_df$fit_method[results_df$fit_method == "gradientTreeBoostCRF"] <- "gradient-tree-CRF-LOP"
results_df$fit_method[results_df$fit_method == "gradientTreeBoostCRF-worst"] <- "gradient-tree-CRF-worst"
results_df$fit_method[results_df$fit_method == "gradientTreeBoostCRF-med"] <- "gradient-tree-CRF-med"
results_df$fit_method[results_df$fit_method == "gradientTreeBoostCRF-best"] <- "gradient-tree-CRF-best"
results_df$fit_method[results_df$fit_method == "unregularizedParametricCRF"] <- "Par-CRF"
results_df$fit_method[results_df$fit_method == "RFHMM"] <- "RF-HMM"
results_df$fit_method[results_df$fit_method == "MLRHMM"] <- "MLR-HMM"
results_df$fit_method[results_df$fit_method == "RF"] <- "RF"
results_df$fit_method[results_df$fit_method == "SVM"] <- "SVM"
results_df$fit_method[results_df$fit_method == "normalHMM"] <- "Gaussian-mixt-HMM"
results_df$fit_method[results_df$fit_method == "BayesRuleTestParams"] <- "Bayes-Rule"
results_df$fit_method <- factor(results_df$fit_method, levels = c("BB-Par-CRF", "BB-Nonpar-CRF", "RF-CRF", "RF-seq-CRF", "gradient-tree-CRF-LOP", "gradient-tree-CRF-worst", "gradient-tree-CRF-med", "gradient-tree-CRF-best", "Par-CRF", "RF-HMM", "MLR-HMM", "RF", "SVM", "Gaussian-mixt-HMM", "Bayes-Rule"))

simStudyResults <- results_df
for(colind in seq_len(5))
	simStudyResults[, colind] <- as.factor(simStudyResults[, colind])

save(simStudyResults, file = paste0(combined_results_dir, "simStudyResults.rdata"))
