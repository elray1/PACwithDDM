rm(list = ls())
library("rayThesisSimStudy")
library("ggplot2")
library("grid")

L_subject_numbers <- c("22")
R_subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "23", "24", "27", "32", "33", "34")
subject_numbers <- c(L_subject_numbers, R_subject_numbers)

locations <- c("ankle", "hip", "wrist")

DO_file_path_base <- "C:/Kinesiology/JefferData/Free Living data/ACE DO data/ACE"
DO_file_path_termination <- "DO.txt"

DO_adjustment_file_path_base <- "C:/Kinesiology/JefferData/Free Living data/ACE DO data/DO_adjustment_subj"
DO_adjustment_file_path_termination <- ".txt"

for(subj in subject_numbers) {
  DO_file_path <- paste0(DO_file_path_base, subj, DO_file_path_termination)
  
  for(location in locations) {
    DO_adjustment_file_path <- paste0(DO_adjustment_file_path_base, subj, "_", location, DO_adjustment_file_path_termination)
    
    cat('"Overall adjustment:","-1.000"\n', file = DO_adjustment_file_path, append = FALSE)
    
    # Read in the DO data
    raw_DO_data <- read.csv(file = file(DO_file_path), skip = 4, stringsAsFactors = FALSE)
    raw_DO_data <- rbind(raw_DO_data[raw_DO_data$Event.Type == "State start", c(1, 2, 5, 6, 7, 8)], raw_DO_data[dim(raw_DO_data)[1], c(1, 2, 5, 6, 7, 8)])
    colnames(raw_DO_data) <- c("Rel.Time", "Subject", "Behavior", "Behavior.Mod1", "Behavior.Mod2", "Event.Type")
    
    # output adjustments to relative times
    cat('"Original Relative Time (seconds)","Adjustment (seconds)"\n', file = DO_adjustment_file_path, append = TRUE)
    for(row_ind in seq_len(nrow(raw_DO_data))) {
      cat(paste0('"', as.character(raw_DO_data$Rel.Time[row_ind]), '",', '"', "0.000", '"\n'), file = DO_adjustment_file_path, append = TRUE)
    }
  }
}
