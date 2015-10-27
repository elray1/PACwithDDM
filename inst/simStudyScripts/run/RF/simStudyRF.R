rm(list=ls())
options(warn = 2, error = recover)

library("lccrf")
library("rayThesis")
library("randomForest")
library("plyr")
library("rstream")

## Preprocess the data

		fit_data <- list(X = rbind.fill.matrix(lapply(dataset, function(comp) comp$X)),
			y = unlist(lapply(dataset, function(comp) comp[[class_var]])),
			subject = unlist(lapply(seq_along(dataset), function(subj) rep(subj, nrow(dataset[[subj]]$X)))))

		N <- length(dataset)
		D <- ncol(fit_data$X)

		if(identical(location, "ankle")) {
			if(identical(class_var, "y_category4")) {
				rng_seed <- 597247
			} else if(identical(class_var, "y_intensity")) {
				rng_seed <- 541479
			}
		} else if(identical(location, "hip")) {
			if(identical(class_var, "y_category4")) {
				rng_seed <- 987785
			} else if(identical(class_var, "y_intensity")) {
				rng_seed <- 598199
			}
		} else if(identical(location, "wrist")) {
			if(identical(class_var, "y_category4")) {
				rng_seed <- 722889
			} else if(identical(class_var, "y_intensity")) {
				rng_seed <- 415681
			}
		}

		# Get rng substream offset corresponding to the combination of fit_method, reduced_trans_mat_parameterization, and update_trans:
		# Build a data frame with the number of substreams used for each combination,
		# including methods that are not run by this code, but that do use rng substreams.
		# The substream offset is then 1 + the total number of substreams used by all earlier methods

		# possible levels for method
		all_fit_methods <- c("RFCRF", "RFCRFseqbag", "parametricBoostCRF", "RF", "normalHMM", "L2RegularizedCRF",
			"gradientTreeBoostCRF", "baggedFeatureSubsetGradientTreeBoostCRF", "MLRHMM", "RFHMM")

		# whether the reduced_trans_mat_parameterization and update_trans parameters are relevant to each method
		reduced_trans_mat_parameterization_relevant <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
		update_trans_relevant <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)

		# create the data frame discussed above -- initialize with all combinations, then remove irrelevant combinations
		substream_df <- data.frame(fit_method = rep(all_fit_methods, each = 4),
			reduced_trans = rep(c(TRUE, FALSE), each = 2, times = length(all_fit_methods)),
			update_trans = rep(c(TRUE, FALSE), times = 2 * length(all_fit_methods)))

		# remove rows where reduced_trans_mat_parameterization is not relevant and reduced_trans is FALSE
		to_remove <- substream_df$fit_method %in% all_fit_methods[!reduced_trans_mat_parameterization_relevant] & !substream_df$reduced_trans
		substream_df <- substream_df[!to_remove, ]

		# remove rows where update_trans is not relevant and update_trans is FALSE
		to_remove <- substream_df$fit_method %in% all_fit_methods[!update_trans_relevant] & !substream_df$update_trans
		substream_df <- substream_df[!to_remove, ]

		# set number of substreams used -- depends on fit_method only
		num_substreams_used <- c(N * 1001, # RFCRF
			N * 1001, # RFCRFseqbag
			N * 1001, # parametricBoostCRF: 1 to generate sequence bag groups + M_bag = 1000, one for each bag group
			1, # RF: 1 stream for all subjects
			2, # normalHMM: 1 stream for all subjects in each case of FullTrans and ReducedTrans
			N * 11, # L2RegularizedCRF: 1 to generate sequence crossvalidation groups + K_crossval = 10, one for each crossvalidation group
			N * (2 + 10 * D * 30 + 1000), # gradientTreeBoostCRF: 1 to generate sequence bag groups, 1 + crf_control$K_crossval * D * 30 to perform cross-validation to select tree depth, M_bag = 1000, one for each bag group
			N * (2 + 10 * D * 30 + 1000), # baggedFeatureSubsetGradientTreeBoostCRF
			0, # MLRHMM: the method is not random
			0) # RFHMM: the method uses the fits created for the RF method above.

		substream_df$substreams_used <- sapply(substream_df$fit_method, function(fm) {
			num_substreams_used[all_fit_methods == fm]
		})

		current_scenario_row_ind <- which(substream_df$fit_method == "RF")[[1]]

		rstream_substream_offset <- 1 + sum(substream_df$substreams_used[seq_len(current_scenario_row_ind - 1)])

		fit_data$y <- factor(fit_data$y)

		combined_results <- vector("list", 33)

		set.seed(rng_seed)
		rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

		# advance to correct substream
		rstream.reset(rngstream)
		for(i in seq_len(rstream_substream_offset))
		  rstream.nextsubstream(rngstream)

		# set rstream as the RNG used by R, using the rngstream object
		rstream.RNG(rngstream)


		for(subj in seq_len(33)) {
			train_subjects <- unique(fit_data$subject)
			train_subjects <- train_subjects[train_subjects != subj]
			T_train <- as.vector(table(fit_data$subject[fit_data$subject != subj]))
  
			train_data <- data.frame(fit_data$X[fit_data$subject != subj, , drop = FALSE])
			train_data$y <- fit_data$y[fit_data$subject != subj]
  
			test_data <- data.frame(fit_data$X[fit_data$subject == subj, , drop = FALSE])
			test_data$y <- fit_data$y[fit_data$subject == subj]
  
			fit <- randomForest(y ~ ., data = train_data)
			y_pred <- predict(fit, test_data[, !(colnames(test_data) %in% c("y", "Y"))])
			log_class_probs <- log(predict(fit, test_data[, !(colnames(test_data) %in% c("y", "Y"))], type = "prob"))
  
			confusion_matrix <- table(test_data$y, y_pred)
  
			num_correct <- sum(y_pred == test_data$y)
			prop_correct <- num_correct / length(test_data$y)
  
			combined_results[[subj]] <- list(subj = subj, fit = fit, y_pred = y_pred, log_class_probs = log_class_probs, confusion_matrix = confusion_matrix, prop_correct = prop_correct)
		}

		sapply(combined_results, function(comp) comp$prop_correct)

		save(combined_results, file = paste0("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_", location, "_", class_var, "_RF.Rdata"))
	}
}
