rm(list = ls())
library("rayThesis")
library("irr")

results_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Mannini/results"
combined_results_dir <- "C:/Stat/HMM/HMMEnsembles/rayThesis/data/"

N <- 33

na_vec <- rep(NA, 2 * 1 * (8 * 2 * 33 + 1 * 1 * 33))
results_df <- data.frame(location = na_vec, class_var = na_vec, fit_method = na_vec, trans_mat_parameterization = na_vec, subject = na_vec,
	prop_correct = na_vec, F1_score_macro = na_vec, mse_pred = na_vec)

#results_df[, 6:7] <- as.numeric(results_df[, 6:7])

#	prop_correct = na_vec, kappa = na_vec, F1_score_micro = na_vec, F1_score_macro = na_vec)
row_num <- 1

for(location in c("ankle", "wrist")) {
	if(identical(location, "ankle")) {
		fit_data <- ManniniAnkleCorrected
	} else if(identical(location, "wrist")) {
		fit_data <- ManniniWristCorrected
	}
	for(class_var in c("y_category4")) {
		for(fit_method in c("baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "RFCRF", "RFCRFseqbag", "gradientTreeBoostCRF", "unregularizedParametricCRF", "RF", "SVM", "RFHMM", "MLRHMM", "normalHMM")) {
			if(fit_method %in% c("RF", "SVM")) {
				trans_mat_parameterization_levels <- "Reduced"
			} else if(fit_method %in% c("unregularizedParametricCRF")) {
				trans_mat_parameterization_levels <- "Full"
			} else {
				trans_mat_parameterization_levels <- c("Full", "Reduced")
			}
			for(trans_mat_parameterization in trans_mat_parameterization_levels) {
				for(subject in seq_len(N)) {
					results_filename <- paste0(trans_mat_parameterization, "Trans_subject", subject, ".Rdata")

					load(file.path(results_dir, location, class_var, fit_method, results_filename))

					num_classes <- 4

					if(identical(fit_method, "gradientTreeBoostCRF")) {
						F1_score_macro_by_m <- lapply(confusion_matrix_by_m, function(confusion_matrix, num_classes) {
							temp <- as.data.frame(confusion_matrix)
							inds <- as.matrix(temp[, 1:2])
							storage.mode(inds) <- "integer"
							vals <- temp[, 3]

							confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
							confusion_matrix[inds] <- temp[, 3]

		#					# calculate kappa (relative to true values)
		#					kappa <- kappa2(cbind(fit_data[[subject]][[class_var]], y_pred))$value
		#
							tp_by_class <- diag(confusion_matrix)
							tn_by_class <- sapply(seq_along(tp_by_class), function(ind) sum(tp_by_class[-ind]))
							fp_by_class <- apply(confusion_matrix, 2, sum) - tp_by_class
							fn_by_class <- apply(confusion_matrix, 1, sum) - tp_by_class
		#
		#					# calculate micro F1
		#					precision_micro <- sum(tp_by_class) / (sum(tp_by_class) + sum(fp_by_class))
		#					recall_micro <- sum(tp_by_class) / (sum(tp_by_class) + sum(fn_by_class))
		#					F1_score_micro <- 2 * precision_micro * recall_micro / (precision_micro + recall_micro)
		#
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
								est_class_probs <- exp(log_class_probs[[1]])
							} else if(fit_method %in% c("SVM")) {
								est_class_probs <- matrix(NA, nrow = length(y_pred), ncol = num_classes)
							} else {
								est_class_probs <- exp(log_class_probs)
							}
							y_true_as_matrix <- matrix(0, nrow = length(fit_data[[subject]][[class_var]]), ncol = num_classes)
							y_true_as_matrix[cbind(seq_along(fit_data[[subject]][[class_var]]), as.integer(fit_data[[subject]][[class_var]]))] <- 1
							mse_pred <- sum((y_true_as_matrix - est_class_probs)^2) / length(fit_data[[subject]][[class_var]])

							return(mse_pred)
						})

						# worst values
						results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = "gradientTreeBoostCRF-worst",
							trans_mat_parameterization = trans_mat_parameterization, subject = subject,
							prop_correct = min(unlist(prop_correct_by_m)), F1_score_macro = min(unlist(F1_score_macro_by_m)), mse_pred = max(unlist(mse_pred_by_m)))
						row_num <- row_num + 1

						# median values
						results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = "gradientTreeBoostCRF-med",
							trans_mat_parameterization = trans_mat_parameterization, subject = subject,
							prop_correct = median(unlist(prop_correct_by_m)), F1_score_macro = median(unlist(F1_score_macro_by_m)), mse_pred = median(unlist(mse_pred_by_m)))
						row_num <- row_num + 1

						# best values
						results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = "gradientTreeBoostCRF-best",
							trans_mat_parameterization = trans_mat_parameterization, subject = subject,
							prop_correct = max(unlist(prop_correct_by_m)), F1_score_macro = max(unlist(F1_score_macro_by_m)), mse_pred = min(unlist(mse_pred_by_m)))
						row_num <- row_num + 1
					}

					temp <- as.data.frame(confusion_matrix)
					inds <- as.matrix(temp[, 1:2])
					storage.mode(inds) <- "integer"
					vals <- temp[, 3]

					confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
					confusion_matrix[inds] <- temp[, 3]

#					# calculate kappa (relative to true values)
#					kappa <- kappa2(cbind(fit_data[[subject]][[class_var]], y_pred))$value
#
					tp_by_class <- diag(confusion_matrix)
					tn_by_class <- sapply(seq_along(tp_by_class), function(ind) sum(tp_by_class[-ind]))
					fp_by_class <- apply(confusion_matrix, 2, sum) - tp_by_class
					fn_by_class <- apply(confusion_matrix, 1, sum) - tp_by_class
#
#					# calculate micro F1
#					precision_micro <- sum(tp_by_class) / (sum(tp_by_class) + sum(fp_by_class))
#					recall_micro <- sum(tp_by_class) / (sum(tp_by_class) + sum(fn_by_class))
#					F1_score_micro <- 2 * precision_micro * recall_micro / (precision_micro + recall_micro)
#
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
						est_class_probs <- exp(log_class_probs[[1]])
					} else if(fit_method %in% c("SVM")) {
						est_class_probs <- matrix(NA, nrow = length(y_pred), ncol = num_classes)
					} else {
						est_class_probs <- exp(log_class_probs)
					}
					y_true_as_matrix <- matrix(0, nrow = length(fit_data[[subject]][[class_var]]), ncol = num_classes)
					y_true_as_matrix[cbind(seq_along(fit_data[[subject]][[class_var]]), as.integer(fit_data[[subject]][[class_var]]))] <- 1
					mse_pred <- sum((y_true_as_matrix - est_class_probs)^2) / length(fit_data[[subject]][[class_var]])

					results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = fit_method,
						trans_mat_parameterization = trans_mat_parameterization, subject = subject,
						prop_correct = prop_correct, F1_score_macro = F1_score_macro, mse_pred = mse_pred)
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
results_df$fit_method <- factor(results_df$fit_method, levels = c("BB-Par-CRF", "BB-Nonpar-CRF", "RF-CRF", "RF-seq-CRF", "gradient-tree-CRF-LOP", "gradient-tree-CRF-worst", "gradient-tree-CRF-med", "gradient-tree-CRF-best", "Par-CRF", "RF-HMM", "MLR-HMM", "RF", "SVM", "Gaussian-mixt-HMM"))

results_df$location <- as.character(results_df$location)
results_df$location[results_df$location == "ankle"] <- "Ankle"
results_df$location[results_df$location == "wrist"] <- "Wrist"

ManniniTypeResults <- results_df
for(colind in seq_len(5))
	ManniniTypeResults[, colind] <- as.factor(ManniniTypeResults[, colind])

save(ManniniTypeResults, file = paste0(combined_results_dir, "ManniniTypeResults.rdata"))
