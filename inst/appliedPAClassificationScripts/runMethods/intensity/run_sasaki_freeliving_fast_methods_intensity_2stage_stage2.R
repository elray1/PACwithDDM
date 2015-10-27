options(warn = 2)


class_var <- "y_intensity"
for(setting in c("freeliving")) {
	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
		for(fit_method in c("RF", "RFHMM")) {
#		for(fit_method in c("RFHMM")) {
			setwd(file.path("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/intensity/results", setting, "2stage", "stage2", location, class_var, fit_method))

			if(fit_method %in% c("RFHMM", "MLRHMM", "normalHMM")) {
#				reduced_trans_mat_parameterization_levels <- c(TRUE, FALSE)
#				reduced_trans_mat_parameterization_levels <- TRUE
				reduced_trans_mat_parameterization_levels <- FALSE
			} else {
				reduced_trans_mat_parameterization_levels <- TRUE
			}
			for(reduced_trans_mat_parameterization in reduced_trans_mat_parameterization_levels) {
				update_trans_levels <- TRUE
				for(update_trans in update_trans_levels) {
					
					if(identical(setting, "freeliving")) {
						N <- 15
					} else if(identical(setting, "lab")) {
						if(identical(location, "ankle")) {
							N <- 34
						} else {
							N <- 35
						}
					}
					
					for(subject in seq(from = 1, to = N)) {
#					for(subject in 2) {
						system_cmd <- paste0("R --vanilla --args ", subject, " ", setting, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization, " ", update_trans,
										 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/runMethods/intensity/sasakiLOSOcrossval_2stage_stage2.R ")
					
						system(system_cmd, intern = TRUE)
					}
				}
			}
		}
	}
}
