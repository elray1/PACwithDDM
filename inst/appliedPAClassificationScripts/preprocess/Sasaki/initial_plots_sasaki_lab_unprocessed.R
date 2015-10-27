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

DO_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Lab Data", "DO", "start and stop record ACE.csv")

# Read in the DO data -- used here to obtain actigraph file path
raw_DO_data <- read.csv(file = file(DO_file_path), stringsAsFactors = FALSE)

demographics_file_path <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/Demographics ACE.csv"

sfExportAll()


#subject_numbers <- subject_numbers[1]

#debug(Sasaki_combine_lab_actigraph_and_DO_files)
sfLapply(subject_numbers, function(subj) {
#lapply(subject_numbers, function(subj) {
#  for(location in c("ankle", "hip", "wrist")) {
  for(location in c("hip")) {
    if(identical(location, "ankle")) {
      location_first_upper <- "Ankle"
    } else if(identical(location, "hip")) {
      location_first_upper <- "Hip"
    } else if(identical(location, "wrist")) {
      location_first_upper <- "Wrist"
    } else {
      stop("Invalid location")
    }

	oxy_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/ACE Oxycon/Oxycon Files/", raw_DO_data$File.name.Oxy[which.max(raw_DO_data$Subject.ID == paste0("ACE", subj) & raw_DO_data$File.name.Oxy != "")])
      
    temp_plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/lab/plots/initial/unprocessed/", location, "/temp_", subj, ".pdf")
    
    plot_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/lab/plots/initial/unprocessed/", location, "/unprocessed_data_classification_plot_behavior_subj", subj, ".png")

    actigraph_file_name <- paste0(raw_DO_data[which(substring(raw_DO_data[, 1], 4, 5) == subj)[1], 7], "RAW.CSV")
    actigraph_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Lab data", "csv", actigraph_file_name)
    
    #			DO_adjustment_file_name <- paste0("DO_adjustment_subj", subj, ".csv")
    #			DO_adjustment_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Free Living data", "rayDOadjustments", DO_adjustment_file_name)
    DO_adjustment_file_path_base <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/DO/DO_adjustment_subj"
    DO_adjustment_file_path_termination <- ".csv"
    DO_adjustment_file_path <- paste0(DO_adjustment_file_path_base, subj, DO_adjustment_file_path_termination)
    
    unprocessed_data <- Sasaki_combine_lab_actigraph_and_DO_files(subj, actigraph_file_path, DO_file_path, DO_adjustment_file_path, oxy_file_path = oxy_file_path, demographics_file_path = demographics_file_path, drop_unlabeled = FALSE)
    unprocessed_data$vm <- apply(as.matrix(unprocessed_data[, c("x", "y", "z")]), 1, function(row_dat) sqrt(sum(row_dat^2)))
    
    unprocessed_data$behavior <- as.character(unprocessed_data$behavior)
    unprocessed_data$behavior[unprocessed_data$behavior == "Simulated Bowling"] <- "Bowling"
    unprocessed_data$behavior[unprocessed_data$behavior == "Slow walk (1.8 mph)"] <- "Slow walk"
    unprocessed_data$behavior[unprocessed_data$behavior == "Crossword puzzles"] <- "Crossword"
    unprocessed_data$behavior[unprocessed_data$behavior == "Playing cards"] <- "Cards"
    unprocessed_data$behavior[unprocessed_data$behavior == "Self-care (miscellaneous)"] <- "Self-care"
    unprocessed_data$behavior[unprocessed_data$behavior == "Organizing the room"] <- "Organizing"
    unprocessed_data$behavior[unprocessed_data$behavior == "Carrying groceries"] <- "Groceries"
    unprocessed_data$behavior <- factor(unprocessed_data$behavior, levels = c("Standing", "Lying down", "Seated", "Crossword", "Cards", "Dusting", "Gardening", "Vacuuming", "Self-care", "Laundry", "Organizing", "Slow walk", "400m walk", "Groceries", "Tai-chi", "Bowling", "Unlabeled"))
    
    plot_main_title <- paste0("Sasaki Lab ", location_first_upper, " Data Subject ", subj)
    
    #debug(make_one_signal_classification_plot)
#    recorded_plots <- make_signal_classification_plots(signal_vars = c("vm", "x", "y", "z"), class_vars = c("behavior", "behavior_mod_2"), data = unprocessed_data,
		recorded_plots <- make_signal_classification_plots(signal_vars = c("vm"), class_vars = c("behavior"), data = unprocessed_data,
																											 sampling_freq = 80, subsample_rate = 2, panel_minutes = 30, main_title = plot_main_title,
                                                       class_var_labels = c("Behavior", "Intensity"), class_var_palette_types = c("manual", "div"),
                                                       class_var_palettes = list(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66FF33", "#FFFF00", "#FFFFFF", "#444444", "#AAAAAA"),
#                                                                                 "RdYlBu"), font_size = 44, line_height = 4, plot_file_path = plot_file_path)
                                                                                 "RdYlBu"), font_size = 44, line_height = 4, panel_legend_relative_width = 5, plot_file_path = plot_file_path, temp_plot_file_path = temp_plot_file_path)
    
    #png(file = plot_file_path, width = 30, height = 20, units = "in", res = 300, type = "cairo-png")
#    pdf(file = plot_file_path, width = 30, height = 20)
#    
#    lapply(recorded_plots, replayPlot)
#    
#    dev.off()
  }
})

sfStop()
