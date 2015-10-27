rm(list = ls())

library("rayThesis")

L_subject_numbers <- c("02", "07", "22", "28")
R_subject_numbers <- c("01", "03", "04", "05", "06", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24", "25", "26", "27", "29", "32", "33", "34", "35", "36", "37")
subject_numbers <- sort(c(L_subject_numbers, R_subject_numbers))
subject_numbers <- "23"

DO_file_path <- file.path(path.package("rayThesis"), "extdata", "Sasaki", "Lab Data", "DO", "start and stop record ACE.csv")

# Read in the DO data
raw_DO_data <- read.csv(file = file(DO_file_path), stringsAsFactors = FALSE)


DO_adjustment_file_path_base <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/DO/DO_adjustment_subj"
DO_adjustment_file_path_termination <- ".csv"

for(subj in subject_numbers) {
  DO_adjustment_file_path <- paste0(DO_adjustment_file_path_base, subj, DO_adjustment_file_path_termination)
  
  cat('"Overall adjustment:","0.000"\n', file = DO_adjustment_file_path, append = FALSE)
  
  DO_data <- raw_DO_data[substring(raw_DO_data[, 1], 4, 5) == subj, ]
  if(identical(subj, "23")) {
	DO_data$Start[1:3] <- c("9:26:00", "9:26:31", "9:27:01")
	DO_data$Stop[1:3] <- c("9:26:30", "9:27:00", "9:27:30")
  }

  
  # output adjustments to times
  cat('"Original Time (seconds)","Adjustment (seconds)"\n', file = DO_adjustment_file_path, append = TRUE)
  for(row_ind in seq_len(nrow(DO_data))) {
  	if(row_ind > 1) {
  		if(!identical(format(as.POSIXlt(strptime(DO_data$Stop[row_ind - 1], "%H:%M:%S")) + 1, "%H:%M:%S"),
  			format(as.POSIXlt(strptime(DO_data$Start[row_ind], "%H:%M:%S")), "%H:%M:%S"))) {
  			cat(paste0('"', as.character(DO_data$Start[row_ind]), '",', '"', "0.000", '"\n'), file = DO_adjustment_file_path, append = TRUE)
  		}
  	} else {
  		cat(paste0('"', as.character(DO_data$Start[row_ind]), '",', '"', "0.000", '"\n'), file = DO_adjustment_file_path, append = TRUE)
  	}
  	cat(paste0('"', as.character(DO_data$Stop[row_ind]), '",', '"', "0.000", '"\n'), file = DO_adjustment_file_path, append = TRUE)
  }
}
