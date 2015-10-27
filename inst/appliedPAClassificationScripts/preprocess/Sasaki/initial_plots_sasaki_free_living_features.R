rm(list = ls())
options(error = recover)

library("snowfall")

sfInit(parallel = TRUE, cpus = 1, type = "SOCK")

sfLibrary("rayThesis", character.only = TRUE)
sfLibrary("ggplot2", character.only = TRUE)
sfLibrary("grid", character.only = TRUE)
sfLibrary("plyr", character.only = TRUE)

rayThesis_location <- find.package("rayThesis")

L_subject_numbers <- c("22")
R_subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "23", "24", "27", "32", "33", "34")
subject_numbers <- c(L_subject_numbers, R_subject_numbers)

sfExportAll()

subject_numbers <- "22"
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
#lapply(subject_numbers, function(subj) {
#  for(location in c("ankle", "hip", "wrist")) {
  for(location in c("hip")) {
    if(identical(location, "ankle")) {
      location_first_upper <- "Ankle"
      temp_plot_data <- SasakiFreeLivingAnkle
    } else if(identical(location, "hip")) {
      location_first_upper <- "Hip"
      temp_plot_data <- SasakiFreeLivingHip
    } else if(identical(location, "wrist")) {
      location_first_upper <- "Wrist"
      temp_plot_data <- SasakiFreeLivingWrist
    } else {
      stop("Invalid location")
    }
    
    temp_plot_data <- temp_plot_data[[which(sapply(temp_plot_data, function(comp) comp$subject == subj))]]
    plot_data <- as.data.frame(temp_plot_data$X)
    plot_data$y_behavior <- temp_plot_data$y_behavior
    plot_data$y_category5 <- temp_plot_data$y_category5
    plot_data$y_category3 <- temp_plot_data$y_category3
    plot_data$y_intensity <- temp_plot_data$y_intensity
    
    plot_feature_names <- colnames(temp_plot_data$X)
    
    for(class_var in c("y_category5", "y_category3", "y_intensity")) {
	    temp_plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/processed/", location, "/temp1_", subj, "_", class_var, ".pdf")
	    plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/processed/", location, "/features_classification_plot_", class_var, "_subj", subj, ".pdf")
	    plot_main_title <- paste0("Sasaki Free Living ", location_first_upper, " Data Subject ", subj)
    
    #debug(make_one_signal_classification_plot)
#    recorded_plots <- make_signal_classification_plots(signal_vars = c("vm", "x", "y", "z"), class_vars = c("behavior", "behavior_mod_2"), data = unprocessed_data,
	    recorded_plots <- make_signal_classification_plots(signal_vars = plot_feature_names, class_vars = class_var, data = plot_data,
	                                              sampling_freq = 1/12.8, subsample_rate = 1, panel_minutes = 30, main_title = plot_main_title,
	                                              class_var_labels = class_var, class_var_palette_types = c("manual", "div"),
	                                              class_var_palettes = list(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF"),
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


for(location in c("hip")) {
	if(identical(location, "ankle")) {
		location_first_upper <- "Ankle"
		temp_plot_data <- SasakiFreeLivingAnkle
	} else if(identical(location, "hip")) {
		location_first_upper <- "Hip"
		temp_plot_data <- SasakiFreeLivingHip
	} else if(identical(location, "wrist")) {
		location_first_upper <- "Wrist"
		temp_plot_data <- SasakiFreeLivingWrist
	} else {
		stop("Invalid location")
	}
	
	plot_data <- as.data.frame(rbind.fill.matrix(lapply(temp_plot_data, function(comp) comp$X)))
	plot_data$y_behavior <- unlist(lapply(temp_plot_data, function(comp) comp$y_behavior))
	plot_data$y_category5 <- unlist(lapply(temp_plot_data, function(comp) comp$y_category5))
	plot_data$y_category3 <- unlist(lapply(temp_plot_data, function(comp) comp$y_category3))
	plot_data$y_intensity <- unlist(lapply(temp_plot_data, function(comp) comp$y_intensity))
	
	plot_feature_names <- colnames(temp_plot_data[[1]]$X)
	
	for(class_var in c("y_category5", "y_category3", "y_intensity")) {
		temp_plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/processed/", location, "/temp1_", class_var, ".pdf")
		plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/freeliving/plots/initial/processed/", location, "/feature_univar_KDEs_classification_plot_", class_var, ".pdf")
		plot_main_title <- paste0("Sasaki Free Living ", location_first_upper, " Data")
	
#    pdf(file = plot_file_path, width = 30, height = 20)
		pdf(file = plot_file_path)

		for(feature_var in plot_feature_names) {
			classes <- unique(plot_data[[class_var]])
			p <- ggplot() +
				geom_density(aes_string(x = feature_var, colour = class_var), data = plot_data) +
				theme_bw()
			print(p)
		}
		
    dev.off()
	}
}
