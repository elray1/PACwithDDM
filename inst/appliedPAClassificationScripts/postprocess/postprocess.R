rm(list = ls())
library("PACwithDDM")
library("irr")

options(warn = 2, error = recover)

pacwithddm_source_location <- "C:/Stat/HMM/PACwithDDM"

#combined_results_dir <- "C:/Stat/HMM/HMMEnsembles/rayThesis/data/"
#results_dir_base <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/"
#combined_results_dir <- "F:/Evan/PACwithDDM-linux/pacwithddm/data"
#results_dir_base <- "F:/Evan/PACwithDDM-linux/pacwithddm/inst/results"
combined_results_dir <- file.path(pacwithddm_source_location, "data")
results_dir_base <- file.path(pacwithddm_source_location, "inst/results")


for(data_set in c("Mannini", "SasakiLab", "SasakiFreeLiving")) {
  location_levels <- c("ankle", "wrist")

  for(class_var_group in c("Type", "Intensity")) {
    results_df <- data.frame(location = NA, class_var = NA, fit_method = NA, case = NA, subject = NA,
      prop_correct = NA, F1_score_macro = NA, mse_pred = NA)
    row_num <- 1
    
    for(location in location_levels) {
      if(identical(data_set, "SasakiFreeLiving")) {
        N <- 15
        
        if(identical(class_var_group, "Type")) {
          class_var_levels <- "y_category3"
        } else {
          class_var_levels <- c("y_intensity")
        }
        
        if(identical(location, "ankle")) {
          fit_data <- SasakiFreeLivingAnkle
        } else if(identical(location, "hip")) {
          fit_data <- SasakiFreeLivingHip
        } else if(identical(location, "wrist")) {
          fit_data <- SasakiFreeLivingWrist
        }
      } else if(identical(data_set, "SasakiLab")) {
        if(identical(location, "ankle")) {
          N <- 34
        } else {
          N <- 35
        }
        
        if(identical(class_var_group, "Type")) {
          class_var_levels <- "y_category3"
        } else {
          class_var_levels <- "y_intensity"
        }
        
        if(identical(location, "ankle")) {
          fit_data <- SasakiLabAnkle
        } else if(identical(location, "hip")) {
          fit_data <- SasakiLabHip
        } else if(identical(location, "wrist")) {
          fit_data <- SasakiLabWrist
        }
      } else if(identical(data_set, "Mannini")) {
        N <- 33
        
        if(identical(class_var_group, "Type")) {
          class_var_levels <- c("y_category4")
        } else {
          class_var_levels <- c("y_intensity")
        }
        
        if(identical(location, "ankle")) {
          fit_data <- ManniniAnkleCorrected
        } else if(identical(location, "wrist")) {
          fit_data <- ManniniWristCorrected
        }
      }
      
      for(class_var in class_var_levels) {
        if(class_var %in% c("y_category3", "y_category4")) {
          num_classes <- 4
        } else if(identical(class_var, "y_intensity")) {
          if(data_set %in% c("SasakiFreeLiving")) {
            num_classes <- 5
          } else {
            num_classes <- 4
          }
        } else {
          num_classes <- 6
        }
        
        all_fit_methods <- c("normalFMM", "normalHMM", "RF", "parametricBoostMLR", "parametricBoostCRF")
        
        results_dir <- file.path(results_dir_base, data_set, location, class_var)
        
        for(fit_method in all_fit_methods) {
          aggregated_confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
          
          combined_y_true_as_matrix <- matrix(NA, nrow = 1, ncol = num_classes)
          combined_est_class_probs <- matrix(NA, nrow = 1, ncol = num_classes)
          for(subject in seq_len(N)) {
            results_filename <- paste0("results_subject", subject, ".Rdata")

            load(file.path(results_dir_base, data_set, location, class_var, fit_method, results_filename))
            
            temp <- as.data.frame(confusion_matrix)
            inds <- as.matrix(temp[, 1:2])
            
#            if(fit_method %in% c("normalFMM", "normalHMM", "parametricBoostMLR", "parametricBoostCRF") && identical(class_var_group, "Type") && data_set %in% c("SasakiLab", "SasakiFreeLiving")) {
#             if(fit_method %in% c("normalFMM", "normalHMM", "parametricBoostMLR", "parametricBoostCRF") && identical(class_var_group, "Intensity") && data_set %in% c("SasakiFreeLiving")) {
#               inds[, 1] <- as.character(as.integer(inds[, 1]) - 1L)
#             } else if(fit_method %in% c("RF") && identical(class_var_group, "Intensity") && data_set %in% c("SasakiFreeLiving")) {
#               inds[, 1] <- as.character(as.integer(inds[, 1]) - 1L)
#               inds[, 2] <- as.character(as.integer(inds[, 2]) - 1L)
#             }
            
            if(identical(data_set, "SasakiLab") && identical(class_var_group, "Intensity")) {
              inds[inds == "5"] <- "4" # "Vigorous" intensity was not a possible level, but this was treated differently in two different variables..........
            }
            
            storage.mode(inds) <- "integer"
            vals <- temp[, 3]
            
            confusion_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
            confusion_matrix[inds] <- temp[, 3]
            aggregated_confusion_matrix <- aggregated_confusion_matrix + confusion_matrix
            
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
            if(fit_method %in% c("parametricBoostCRF", "parametricBoostMLR")) {
              est_class_probs <- exp(log_class_probs[[1]])
            } else {
              est_class_probs <- exp(log_class_probs)
            }
            
            y_true_as_matrix <- matrix(0, nrow = length(fit_data[[subject]][[class_var]]), ncol = num_classes)
            y_true_cols <- as.integer(fit_data[[subject]][[class_var]])
            if(identical(class_var_group, "Intensity") && data_set %in% c("SasakiFreeLiving")) {
              y_true_cols <- y_true_cols - 1L
            }
            y_true_as_matrix[cbind(seq_along(fit_data[[subject]][[class_var]]), y_true_cols)] <- 1
            mse_pred <- sum((y_true_as_matrix - est_class_probs)^2) / length(fit_data[[subject]][[class_var]])
            
            combined_y_true_as_matrix <- rbind(combined_y_true_as_matrix, y_true_as_matrix)
            combined_est_class_probs <- rbind(combined_est_class_probs, est_class_probs)
            
            results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = fit_method,
              case = NA, subject = as.character(subject),
              prop_correct = prop_correct, F1_score_macro = F1_score_macro, mse_pred = mse_pred)
            row_num <- row_num + 1
          }
          
          ## add summary statistics for aggregated values across subjects
          tp_by_class <- diag(aggregated_confusion_matrix)
          tn_by_class <- sapply(seq_along(tp_by_class), function(ind) sum(tp_by_class[-ind]))
          fp_by_class <- apply(aggregated_confusion_matrix, 2, sum) - tp_by_class
          fn_by_class <- apply(aggregated_confusion_matrix, 1, sum) - tp_by_class
          
          # calculate macro F1
          to_keep <- (apply(aggregated_confusion_matrix, 1, sum) > 0)
          precision_by_class <- tp_by_class[to_keep] / (tp_by_class[to_keep] + fp_by_class[to_keep])
          precision_by_class[is.na(precision_by_class)] <- 0
          recall_by_class <- tp_by_class[to_keep] / (tp_by_class[to_keep] + fn_by_class[to_keep])
          recall_by_class[is.na(recall_by_class)] <- 0
          
          precision_macro <- mean(precision_by_class)
          recall_macro <- mean(recall_by_class)
          F1_score_macro <- 2 * precision_macro * recall_macro / (precision_macro + recall_macro)
          
          # calculate MSE / Brier score
          if(fit_method %in% c("parametricBoostCRF", "parametricBoostMLR")) {
            est_class_probs <- exp(log_class_probs[[1]])
          } else {
            est_class_probs <- exp(log_class_probs)
          }
          if(identical(fit_method, "SVM")) {
            mse_pred <- NA
          } else {
            combined_y_true_as_matrix <- combined_y_true_as_matrix[-1, , drop = FALSE]
            combined_est_class_probs <- combined_est_class_probs[-1, , drop = FALSE]
            mse_pred <- sum((combined_y_true_as_matrix - combined_est_class_probs)^2) / nrow(combined_y_true_as_matrix)
          }
          
          results_df[row_num, ] <- list(location = location, class_var = class_var, fit_method = fit_method,
            case = NA, subject = "Aggregated",
            prop_correct = prop_correct, F1_score_macro = F1_score_macro, mse_pred = mse_pred)
          row_num <- row_num + 1
        }
      }
    }
    
    results_df$fit_method[results_df$fit_method == "normalFMM"] <- "FMM"
    results_df$fit_method[results_df$fit_method == "normalHMM"] <- "HMM"
    results_df$fit_method[results_df$fit_method == "RF"] <- "RF"
    results_df$fit_method[results_df$fit_method == "parametricBoostMLR"] <- "MLR"
    results_df$fit_method[results_df$fit_method == "parametricBoostCRF"] <- "CRF"
    results_df$fit_method <- factor(results_df$fit_method, levels = c("FMM", "HMM", "RF", "MLR", "CRF"))
    
    results_df$location <- as.character(results_df$location)
    results_df$location[results_df$location == "ankle"] <- "Ankle"
    results_df$location[results_df$location == "wrist"] <- "Wrist"
    
    for(colind in seq_len(5))
      results_df[, colind] <- as.factor(results_df[, colind])
    
    save_name <- paste0(data_set, class_var_group, "Results")
    
    assign(save_name, results_df)
    
    save(list = save_name, file = file.path(combined_results_dir, paste0(save_name, ".rdata")))
  }
}
