rm(list = ls())
library("snowfall")

sfInit(parallel = TRUE, cpus = 1, type = "SOCK")

sfLibrary("rayThesis", character.only = TRUE)
sfLibrary("ggplot2", character.only = TRUE)
sfLibrary("grid", character.only = TRUE)

rayThesis_location <- find.package("rayThesis")

L_subject_numbers <- c("22")
R_subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "23", "24", "27", "32", "33", "34")
subject_numbers <- c(L_subject_numbers, R_subject_numbers)

sfExportAll()

#subject_numbers <- "22"
#subject_numbers <- "01"
#subject_numbers <- "04"
#subject_numbers <- "06"
#subject_numbers <- "08"
#subject_numbers <- "11"
#subject_numbers <- "19"
#subject_numbers <- "20"
#subject_numbers <- "21"
#subject_numbers <- subject_numbers[!(subject_numbers %in% c("22", "01"))]

sfLapply(subject_numbers, function(subj) {
#  for(location in c("ankle", "hip", "wrist")) {
  for(location in c("ankle", "wrist")) {
    if(identical(location, "ankle")) {
      location_first_upper <- "Ankle"
    } else if(identical(location, "hip")) {
      location_first_upper <- "Hip"
    } else if(identical(location, "wrist")) {
      location_first_upper <- "Wrist"
    } else {
      stop("Invalid location")
    }
    
    temp_plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/unprocessed/", location, "/temp1_", subj, ".pdf")
    
    plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/unprocessed/", location, "/unprocessed_data_classification_plot_adjusted_DO_behavior_subj", subj, ".png")
    
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
    
    plot_main_title <- paste0("Sasaki Free Living ", location_first_upper, " Data Subject ", subj, " -- DO Times Adjusted")
    
    #debug(make_one_signal_classification_plot)
#    recorded_plots <- make_signal_classification_plots(signal_vars = c("vm", "x", "y", "z"), class_vars = c("behavior", "behavior_mod_2"), data = unprocessed_data,
    recorded_plots <- make_signal_classification_plots(signal_vars = c("vm"), class_vars = c("behavior"), data = unprocessed_data,
                                              sampling_freq = 80, subsample_rate = 1, panel_minutes = 30, main_title = plot_main_title,
                                              class_var_labels = c("Behavior", "Intensity"), class_var_palette_types = c("manual", "div"),
                                              class_var_palettes = list(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF"),
                                                                        "RdYlBu"), font_size = 44, line_height = 4, plot_file_path = plot_file_path, temp_plot_file_path = temp_plot_file_path)
    
    #png(file = plot_file_path, width = 30, height = 20, units = "in", res = 300, type = "cairo-png")
#    pdf(file = plot_file_path, width = 30, height = 20)
#    
#    lapply(recorded_plots, replayPlot)
#    
#    dev.off()
  }
})

sfStop()
