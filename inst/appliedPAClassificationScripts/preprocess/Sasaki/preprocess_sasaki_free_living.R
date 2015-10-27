rm(list = ls())
options(error = recover)

library("snowfall")

sfInit(parallel = TRUE, cpus = 5, type = "SOCK")

sfLibrary("plyr", character.only = TRUE)
sfLibrary("rayThesis", character.only = TRUE)

rayThesis_location <- find.package("rayThesis")

L_subject_numbers <- c("22")
R_subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "23", "24", "27", "32", "33", "34")
subject_numbers <- sort(c(L_subject_numbers, R_subject_numbers))
#subject_numbers <- subject_numbers[5]

sfExportAll()

#debug(col_fft_features)
#debug(Sasaki_preprocess_one_free_living_file)

#for(location in c("ankle", "hip", "wrist")) {
for(location in c("ankle", "wrist")) {
	data_split <- 
	sfLapply(subject_numbers, function(subj, location) {
#	lapply(subject_numbers, function(subj, location) {
#	subj <- subject_numbers[[1]]
		if(identical(location, "ankle")) {
			location_first_upper <- "Ankle"
		} else if(identical(location, "hip")) {
			location_first_upper <- "Hip"
		} else if(identical(location, "wrist")) {
			location_first_upper <- "Wrist"
		} else {
			stop("Invalid location")
		}
		
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
		
		temp <- Sasaki_preprocess_one_free_living_file(actigraph_file_path = actigraph_file_path,
																		DO_file_path = DO_file_path, DO_adjustment_file_path = DO_adjustment_file_path,
																		sampling_freq = 80, window_length = 12.8)
		
		temp$subject <- subj
		temp$location <- location
		
		return(temp)
	}, location)


	if(identical(location, "ankle")) {
		SasakiFreeLivingAnkle <- data_split
		save(SasakiFreeLivingAnkle, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiFreeLivingAnkle.rdata")
	} else if(identical(location, "hip")) {
		SasakiFreeLivingHip <- data_split
		save(SasakiFreeLivingHip, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiFreeLivingHip.rdata")
	} else if(identical(location, "wrist")) {
		SasakiFreeLivingWrist <- data_split
		save(SasakiFreeLivingWrist, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/SasakiFreeLivingWrist.rdata")
	}
}

sfStop()
