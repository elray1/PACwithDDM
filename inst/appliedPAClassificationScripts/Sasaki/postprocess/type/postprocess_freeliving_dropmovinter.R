rm(list = ls())
library("rayThesis")
library("irr")

results_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/results/freeliving_dropMovInter"
combined_results_dir <- "C:/Stat/HMM/HMMEnsembles/rayThesis/data/"

setting <- "freeliving"
N <- 15

results_df <- data.frame(location = NA, class_var = NA, fit_method = NA, trans_mat_parameterization = NA, subject = NA,
	prop_correct = NA, F1_score_macro = NA, mse_pred = NA)


#	prop_correct = na_vec, kappa = na_vec, F1_score_micro = na_vec, F1_score_macro = na_vec)
row_num <- 1

for(location in c("ankle", "hip", "wrist")) {
	if(identical(location, "ankle")) {
		if(identical(setting, "freeliving")) {
			fit_data <- SasakiFreeLivingAnkle
		} else if(identical(setting, "lab")) {
			fit_data <- SasakiLabAnkle
		} else {
			stop("Invalid setting")
		}
	} else if(identical(location, "hip")) {
		if(identical(setting, "freeliving")) {
			fit_data <- SasakiFreeLivingHip
		} else if(identical(setting, "lab")) {
			fit_data <- SasakiLabHip
		} else {
			stop("Invalid setting")
		}
	} else if(identical(location, "wrist")) {
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
#	for(class_var in c("y_category3", "y_category5")) {
	for(class_var in c("y_category5")) {
		for(fit_method in c("baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "RFCRF", "RFCRFseqbag", "gradientTreeBoostCRF", "unregularizedParametricCRF", "RF", "RFHMM", "MLRHMM", "normalHMM")) {
			if(identical(fit_method, "RF")) {
				trans_mat_parameterization_levels <- "Reduced"
			} else if(fit_method %in% c("unregularizedParametricCRF")) {
				trans_mat_parameterization_levels <- "Full"
			} else {
#				trans_mat_parameterization_levels <- c("Full", "Reduced")
				trans_mat_parameterization_levels <- "Full"
			}
			for(trans_mat_parameterization in trans_mat_parameterization_levels) {
				for(subject in seq_len(N)) {
					results_filename <- paste0(trans_mat_parameterization, "Trans_subject", subject, "_dropMovInter.Rdata")

					load(file.path(results_dir, location, class_var, fit_method, results_filename))

					if(identical(class_var, "y_category3")) {
						num_classes <- 4 - 1
					} else if(identical(class_var, "y_category5")) {
						num_classes <- 6 - 1
					}

#					if(identical(fit_method, "gradientTreeBoostCRF")) {
#						selected_fit_ind <- which(unlist(prop_correct_by_m) == median(unlist(prop_correct_by_m)))[1]
#						prop_correct <- median(unlist(prop_correct_by_m))
#						if(is.na(selected_fit_ind)) {
#							selected_fit_ind
#						}
#						confusion_matrix <- confusion_matrix_by_m[selected_fit_ind]
#						y_pred <- y_pred_by_m[selected_fit_ind]
#					}

					temp <- as.data.frame(confusion_matrix)
					inds <- as.matrix(temp[, 1:2])
					storage.mode(inds) <- "integer"
					inds[inds > 3] <- inds[inds > 3] - 1
					vals <- temp[, 3]

					confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
					confusion_matrix[inds] <- temp[, 3]
					confusion_matrix <- confusion_matrix[-3, -3]

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
					y_true_dropMovInter <- fit_data[[subject]][[class_var]]
					y_true_dropMovInter <- y_true_dropMovInter[y_true_dropMovInter != "MovInter"]
					y_true_dropMovInter <-  factor(as.numeric(y_true_dropMovInter), levels = sort(unique(as.numeric(y_true_dropMovInter))))
					y_true_as_matrix <- matrix(0, nrow = length(y_true_dropMovInter), ncol = num_classes)
					y_true_as_matrix[cbind(seq_along(y_true_dropMovInter), as.integer(y_true_dropMovInter))] <- 1
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
results_df$fit_method[results_df$fit_method == "gradientTreeBoostCRF"] <- "gradient-tree-CRF"
results_df$fit_method[results_df$fit_method == "unregularizedParametricCRF"] <- "Par-CRF"
results_df$fit_method[results_df$fit_method == "RFHMM"] <- "RF-HMM"
results_df$fit_method[results_df$fit_method == "MLRHMM"] <- "MLR-HMM"
results_df$fit_method[results_df$fit_method == "RF"] <- "RF"
results_df$fit_method[results_df$fit_method == "normalHMM"] <- "Gaussian-HMM"
results_df$fit_method <- factor(results_df$fit_method, levels = c("BB-Par-CRF", "BB-Nonpar-CRF", "RF-CRF", "RF-seq-CRF", "gradient-tree-CRF", "Par-CRF", "RF-HMM", "MLR-HMM", "RF", "Gaussian-HMM"))

results_df$location <- as.character(results_df$location)
results_df$location[results_df$location == "ankle"] <- "Ankle"
results_df$location[results_df$location == "hip"] <- "Hip"
results_df$location[results_df$location == "wrist"] <- "Wrist"

SasakiFreeLivingTypeDropMovInterResults <- results_df
for(colind in seq_len(5))
	SasakiFreeLivingTypeDropMovInterResults[, colind] <- as.factor(SasakiFreeLivingTypeDropMovInterResults[, colind])

save(SasakiFreeLivingTypeDropMovInterResults, file = paste0(combined_results_dir, "SasakiFreeLivingTypeDropMovInterResults.rdata"))
