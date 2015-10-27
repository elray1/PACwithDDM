rm(list = ls())
library("snowfall")
options(error = recover)
sfInit(parallel = TRUE, cpus = 5, type = "SOCK")

sfLibrary("rayThesis", character.only = TRUE)
sfLibrary("plyr", character.only = TRUE)

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


for(location in c("hip", "wrist", "ankle")) {
	if(identical(location, "ankle")) {
		subject_numbers <- subject_numbers[subject_numbers != "34"]
	}
#for(location in c("hip")) {
	data_split <- 
		sfLapply(subject_numbers, function(subj, location) {
#		lapply(subject_numbers, function(subj, location) {
			#	subj <- subject_numbers[[1]]
			if(identical(location, "ankle")) {
				location_first_upper <- "Ankle"
				actigraph_file_col_num <- 8
			} else if(identical(location, "hip")) {
				location_first_upper <- "Hip"
				actigraph_file_col_num <- 7
			} else if(identical(location, "wrist")) {
				location_first_upper <- "Wrist"
				actigraph_file_col_num <- 6
			} else {
				stop("Invalid location")
			}

			oxy_file_path <- paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/ACE Oxycon/Oxycon Files/", raw_DO_data$File.name.Oxy[which.max(raw_DO_data$Subject.ID == paste0("ACE", subj) & raw_DO_data$File.name.Oxy != "")])
			

			actigraph_file_name <- paste0(raw_DO_data[which(substring(raw_DO_data[, 1], 4, 5) == subj)[1], actigraph_file_col_num], "RAW.CSV")
			actigraph_file_path <- file.path(rayThesis_location, "extdata", "Sasaki", "Lab data", location, "csv", actigraph_file_name)
			
			DO_adjustment_file_path_base <- "C:/Stat/HMM/HMMEnsembles/rayThesis/inst/extdata/Sasaki/Lab Data/DO/DO_adjustment_subj"
			DO_adjustment_file_path_termination <- ".csv"
			DO_adjustment_file_path <- paste0(DO_adjustment_file_path_base, subj, DO_adjustment_file_path_termination)
			
		    temp <- Sasaki_preprocess_one_lab_file(subj = subj, actigraph_file_path = actigraph_file_path,
				DO_file_path = DO_file_path, DO_adjustment_file_path = DO_adjustment_file_path,
				oxy_file_path = oxy_file_path, demographics_file_path = demographics_file_path,
				sampling_freq = 80, window_length = 12.8, drop_unlabeled = TRUE)
			
			temp$subject <- subj
			temp$location <- location
			
			return(temp)
		}, location)
	
	
	if(identical(location, "ankle")) {
		SasakiLabAnkle <- data_split
		save(SasakiLabAnkle, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiLabAnkle.rdata")
	} else if(identical(location, "hip")) {
		SasakiLabHip <- data_split
		save(SasakiLabHip, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiLabHip.rdata")
	} else if(identical(location, "wrist")) {
		SasakiLabWrist <- data_split
		save(SasakiLabWrist, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiLabWrist.rdata")
	}
}

sfStop()
