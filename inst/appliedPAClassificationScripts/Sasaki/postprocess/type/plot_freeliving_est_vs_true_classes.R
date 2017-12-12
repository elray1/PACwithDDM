rm(list = ls())

library("snowfall")

sfInit(parallel = TRUE, cpus = 8, type = "SOCK")

sfLibrary("rayThesis", character.only = TRUE)
sfLibrary("ggplot2", character.only = TRUE)
sfLibrary("grid", character.only = TRUE)
sfLibrary("plyr", character.only = TRUE)

rayThesis_location <- find.package("rayThesis")


fl_L_subject_numbers <- c("22")
fl_R_subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "23", "24", "27", "32", "33", "34")
fl_subject_numbers <- sort(c(fl_L_subject_numbers, fl_R_subject_numbers))

lab_L_subject_numbers <- c("02", "07", "22", "28")
lab_R_subject_numbers <- c("01", "03", "04", "05", "06", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24", "25", "26", "27", "29", "32", "33", "34", "35", "36", "37")
lab_subject_numbers <- sort(c(lab_L_subject_numbers, lab_R_subject_numbers))

sfExportAll()

for(setting in c("freeliving", "lab")) {
	if(identical(setting, "freeliving")) {
		results_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/results/freeliving"
		results_plots_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/estvslabeled"
		L_subject_numbers <- fl_L_subject_numbers
		R_subject_numbers <- fl_R_subject_numbers
		subject_numbers <- fl_subject_numbers
	} else {
		results_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/results/lab"
		results_plots_dir <- "C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/lab/plots/estvslabeled"
		L_subject_numbers <- lab_L_subject_numbers
		R_subject_numbers <- lab_R_subject_numbers
		subject_numbers <- lab_subject_numbers
	}

	junk <- sfLapply(subject_numbers, function(subj, setting, L_subject_numbers, R_subject_numbers, subject_numbers, results_dir, results_plots_dir) {
		for(location in c("ankle", "hip", "wrist")) {
			if(identical(location, "ankle")) {
				location_first_upper <- "Ankle"
				if(identical(setting, "freeliving")) {
					fit_data <- SasakiFreeLivingAnkle
					N <- 15
				} else if(identical(setting, "lab")) {
					fit_data <- SasakiLabAnkle
					actigraph_file_col_num <- 8
					N <- 34
				} else {
					stop("Invalid setting")
				}
			} else if(identical(location, "hip")) {
				location_first_upper <- "Hip"
				if(identical(setting, "freeliving")) {
					fit_data <- SasakiFreeLivingHip
					N <- 15
				} else if(identical(setting, "lab")) {
					fit_data <- SasakiLabHip
					actigraph_file_col_num <- 7
					N <- 35
				} else {
					stop("Invalid setting")
				}
			} else if(identical(location, "wrist")) {
				location_first_upper <- "Wrist"
				if(identical(setting, "freeliving")) {
					fit_data <- SasakiFreeLivingWrist
					N <- 15
				} else if(identical(setting, "lab")) {
					fit_data <- SasakiLabWrist
					actigraph_file_col_num <- 6
					N <- 35
				} else {
					stop("Invalid setting")
				}
			} else {
				stop("Invalid location")
			}
	
	
    
			if(identical(setting, "freeliving")) {
				L_actigraph_file_path_termination <- paste0(substr(location_first_upper, 1, 1), "L80DORAW.csv")
				R_actigraph_file_path_termination <- paste0(substr(location_first_upper, 1, 1), "R80DORAW.csv")
				if(subj %in% L_subject_numbers) {
				  actigraph_file_name <- paste0("AG", subj, L_actigraph_file_path_termination)
				} else {
				  actigraph_file_name <- paste0("AG", subj, R_actigraph_file_path_termination)
				}
				actigraph_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Free Living data", "Actigraph", location_first_upper, "csv", actigraph_file_name)
    
				DO_file_name <- paste0("ACE", subj, "DO.txt")
				DO_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Free Living data", "ACE DO data", DO_file_name)
    
				DO_adjustment_file_name <- paste0("DO_adjustment_subj", subj, "_", "hip", ".txt")
				DO_adjustment_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Free Living data", "rayDOadjustments", DO_adjustment_file_name)
    
				unprocessed_data <- Sasaki_combine_free_living_actigraph_and_DO_files(actigraph_file_path, DO_file_path, DO_adjustment_file_path)
				unprocessed_data$vm <- apply(as.matrix(unprocessed_data[, c("x", "y", "z")]), 1, function(row_dat) sqrt(sum(row_dat^2)))

				labeled_summarized_data <- Sasaki_preprocess_one_free_living_file(actigraph_file_path = actigraph_file_path,
																	DO_file_path = DO_file_path, DO_adjustment_file_path = DO_adjustment_file_path,
																	sampling_freq = 80, window_length = 12.8)
			} else {
				demographics_file_path <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/Demographics ACE.csv"

				DO_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Lab Data", "DO", "start and stop record ACE.csv")

				# Read in the DO data -- used here to obtain actigraph file path
				raw_DO_data <- read.csv(file = file(DO_file_path), stringsAsFactors = FALSE)

				oxy_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/ACE Oxycon/Oxycon Files/", raw_DO_data$File.name.Oxy[which.max(raw_DO_data$Subject.ID == paste0("ACE", subj) & raw_DO_data$File.name.Oxy != "")])
			

				actigraph_file_name <- paste0(raw_DO_data[which(substring(raw_DO_data[, 1], 4, 5) == subj)[1], actigraph_file_col_num], "RAW.CSV")
				actigraph_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Lab data", location, "csv", actigraph_file_name)
			
				DO_adjustment_file_path_base <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/DO/DO_adjustment_subj"
				DO_adjustment_file_path_termination <- ".csv"
				DO_adjustment_file_path <- paste0(DO_adjustment_file_path_base, subj, DO_adjustment_file_path_termination)

			    unprocessed_data <- Sasaki_combine_lab_actigraph_and_DO_files(subj, actigraph_file_path, DO_file_path, DO_adjustment_file_path, oxy_file_path = oxy_file_path, demographics_file_path = demographics_file_path, drop_unlabeled = TRUE)
			    unprocessed_data$vm <- apply(as.matrix(unprocessed_data[, c("x", "y", "z")]), 1, function(row_dat) sqrt(sum(row_dat^2)))

				labeled_summarized_data <- Sasaki_preprocess_one_lab_file(subj = subj, actigraph_file_path = actigraph_file_path,
					DO_file_path = DO_file_path, DO_adjustment_file_path = DO_adjustment_file_path,
					oxy_file_path = oxy_file_path, demographics_file_path = demographics_file_path,
					sampling_freq = 80, window_length = 12.8, drop_unlabeled = TRUE)
			}
		
			unprocessed_data$y_category3_labeled <- rep(labeled_summarized_data$y_category3, each = 80 * 12.8)[seq_along(unprocessed_data$vm)]
			unprocessed_data$y_category5_labeled <- rep(labeled_summarized_data$y_category5, each = 80 * 12.8)[seq_along(unprocessed_data$vm)]


#			for(fit_method in c("RFCRF", "RFCRFseqbag", "baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "gradientTreeBoostCRF", "RF", "normalHMM", "RFHMM", "MLRHMM")) {
			for(fit_method in c("RFCRF")) {
				if(identical(fit_method, "RF")) {
					trans_mat_parameterization_levels <- "Reduced"
				} else {
					trans_mat_parameterization_levels <- "Full"
				}
	
				results_loadenv <- new.env()
				results_filename <- paste0("FullTrans_subject", which(subject_numbers == subj), ".Rdata")
				load(file.path(results_dir, location, "y_category3", fit_method, results_filename), envir = results_loadenv)
				unprocessed_data$y_category3_estimated <- factor(rep(levels(labeled_summarized_data$y_category3)[results_loadenv$y_pred], each = 80 * 12.8)[seq_along(unprocessed_data$vm)], levels = levels(labeled_summarized_data$y_category3))

				load(file.path(results_dir, location, "y_category5", fit_method, results_filename), envir = results_loadenv)
				unprocessed_data$y_category5_estimated <- factor(rep(levels(labeled_summarized_data$y_category5)[results_loadenv$y_pred], each = 80 * 12.8)[seq_along(unprocessed_data$vm)], levels = levels(labeled_summarized_data$y_category5))

				for(class_var in c("y_category3", "y_category5")) {
					if(identical(setting, "freeliving")) {
						if(identical(class_var, "y_category3")) {
							plot_main_title <- paste0("Sasaki Free Living ", location_first_upper, " Data Subject ", which(subject_numbers == subj), ", 4 Classes\nLabeled (Bottom) vs Estimated (Top) Class, Method: ", fit_method)
							temp_plot_file_name <- paste0("tempplot_4classes_", which(subject_numbers == subj), ".png")
							plot_file_name <- paste0("plot_4classes_", which(subject_numbers == subj), ".png")
						} else {
							plot_main_title <- paste0("Sasaki Free Living ", location_first_upper, " Data Subject ", which(subject_numbers == subj), ", 6 Classes\nLabeled (Bottom) vs Estimated (Top) Class, Method: ", fit_method)
							temp_plot_file_name <- paste0("tempplot_6classes_", which(subject_numbers == subj), ".png")
							plot_file_name <- paste0("plot_6classes_", which(subject_numbers == subj), ".png")
						}
					} else {
						if(identical(class_var, "y_category3")) {
							plot_main_title <- paste0("Sasaki Lab ", location_first_upper, " Data Subject ", which(subject_numbers == subj), ", 4 Classes\nLabeled (Bottom) vs Estimated (Top) Class, Method: ", fit_method)
							temp_plot_file_name <- paste0("tempplot_4classes_", which(subject_numbers == subj), ".png")
							plot_file_name <- paste0("plot_4classes_", which(subject_numbers == subj), ".png")
						} else {
							plot_main_title <- paste0("Sasaki Lab ", location_first_upper, " Data Subject ", which(subject_numbers == subj), ", 6 Classes\nLabeled (Bottom) vs Estimated (Top) Class, Method: ", fit_method)
							temp_plot_file_name <- paste0("tempplot_6classes_", which(subject_numbers == subj), ".png")
							plot_file_name <- paste0("plot_6classes_", which(subject_numbers == subj), ".png")
						}
					}

					temp_plot_file_path <- file.path(results_plots_dir, location, class_var, fit_method, temp_plot_file_name)
					plot_file_path <- file.path(results_plots_dir, location, class_var, fit_method, plot_file_name)
    

#					undebug(make_one_signal_multi_classification_plot)
				    recorded_plots <- make_signal_multi_classification_plots(signal_vars = c("vm"), class_vars = paste0(class_var, c("_labeled", "_estimated")), data = unprocessed_data,
						sampling_freq = 80, subsample_rate = 1, panel_minutes = 30, main_title = plot_main_title,
						class_var_labels = c("Behavior"), class_var_palette_types = c("manual"),
						class_var_palettes = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF"), font_size = 36, line_height = 8, panel_legend_relative_width = 4, plot_file_path = plot_file_path, temp_plot_file_path = temp_plot_file_path)
#						class_var_palettes = list(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF"),
#							"RdYlBu"), font_size = 44, line_height = 4, plot_file_path = plot_file_path, temp_plot_file_path = temp_plot_file_path)
				}
			}
		}
	}, setting = setting, L_subject_numbers = L_subject_numbers, R_subject_numbers = R_subject_numbers, subject_numbers = subject_numbers, results_dir = results_dir, results_plots_dir = results_plots_dir)
}
