rm(list = ls())
library("snowfall")
options(error = recover)
sfInit(parallel = TRUE, cpus = 8, type = "SOCK")

sfLibrary("rayThesis", character.only = TRUE)
sfLibrary("ggplot2", character.only = TRUE)
sfLibrary("grid", character.only = TRUE)

L_subject_numbers <- c("02", "07", "22", "28")
R_subject_numbers <- c("01", "03", "04", "05", "06", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24", "25", "26", "27", "29", "32", "33", "34", "35", "36", "37")
subject_numbers <- sort(c(L_subject_numbers, R_subject_numbers))
#subject_numbers <- subject_numbers[3]
#subject_numbers <- "26"

rayThesis_location <- path.package("rayThesis")

sfExportAll()


#subject_numbers <- subject_numbers[1]

#debug(Sasaki_combine_lab_actigraph_and_DO_files)
sfLapply(subject_numbers, function(subj) {
#lapply(subject_numbers, function(subj) {
#  for(location in c("hip", "wrist", "ankle")) {
  for(location in c("wrist", "ankle")) {
    if(identical(location, "ankle")) {
      subject_numbers <- subject_numbers[subject_numbers != "34"]
      location_first_upper <- "Ankle"
      temp_plot_data <- SasakiLabAnkle
    } else if(identical(location, "hip")) {
      location_first_upper <- "Hip"
      temp_plot_data <- SasakiLabHip
    } else if(identical(location, "wrist")) {
      location_first_upper <- "Wrist"
      temp_plot_data <- SasakiLabWrist
    } else {
      stop("Invalid location")
    }
    
    temp_plot_data <- temp_plot_data[[which(sapply(temp_plot_data, function(comp) comp$subject == subj))]]
    plot_data <- as.data.frame(temp_plot_data$X)
    plot_data$y_behavior <- temp_plot_data$y_behavior
    plot_data$y_category5 <- temp_plot_data$y_category5
#    plot_data$y_category3 <- temp_plot_data$y_category3
    plot_data$y_intensity <- temp_plot_data$y_intensity
      
    plot_feature_names <- colnames(temp_plot_data$X)

    for(class_var in c("y_behavior", "y_category5", "y_intensity")) {
	    temp_plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/lab/plots/initial/processed/", location, "/temp1_", subj, "_", class_var, ".pdf")
	    plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/lab/plots/initial/processed/", location, "/features_classification_plot_", class_var, "_subj", subj, ".pdf")
	    plot_main_title <- paste0("Sasaki Lab ", location_first_upper, " Data Subject ", subj)
    
    #debug(make_one_signal_classification_plot)
#    recorded_plots <- make_signal_classification_plots(signal_vars = c("vm", "x", "y", "z"), class_vars = c("behavior", "behavior_mod_2"), data = unprocessed_data,
	    recorded_plots <- make_signal_classification_plots(signal_vars = plot_feature_names, class_vars = class_var, data = plot_data,
	                                              sampling_freq = 1/12.8, subsample_rate = 1, panel_minutes = 30, main_title = plot_main_title,
	                                              class_var_labels = class_var, class_var_palette_types = c("manual", "div"),
	                                              class_var_palettes = list(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF", "#444444", "#AAAAAA"),
	                                                                        "RdYlBu"), font_size = 44, line_height = 4, plot_file_path = plot_file_path, temp_plot_file_path = temp_plot_file_path)
    
    #png(file = plot_file_path, width = 30, height = 20, units = "in", res = 300, type = "cairo-png")
#    pdf(file = plot_file_path, width = 30, height = 20)
#    
#    lapply(recorded_plots, replayPlot)
#    
#    dev.off()
		}
  }
})


sfStop()
