for(setting in c("lab")) {
	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
		for(class_var in c("y_category3", "y_category5")) {
#		for(class_var in c("y_category3")) {
#			for(fit_method in c("unregularizedParametricCRF", "RF", "RFHMM", "MLRHMM", "normalHMM")) {
			for(fit_method in c("normalHMM")) {
				setwd(file.path("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/results", setting, location, class_var, fit_method))

				if(fit_method %in% c("RFHMM", "MLRHMM", "normalHMM")) {
#					reduced_trans_mat_parameterization_levels <- c(TRUE, FALSE)
					reduced_trans_mat_parameterization_levels <- FALSE
				} else if(fit_method %in% c("unregularizedParametricCRF")) {
					reduced_trans_mat_parameterization_levels <- FALSE
				} else {
					reduced_trans_mat_parameterization_levels <- TRUE
				}
				for(reduced_trans_mat_parameterization in reduced_trans_mat_parameterization_levels) {
					update_trans_levels <- TRUE
					for(update_trans in update_trans_levels) {
						
						if(identical(setting, "freeliving")) {
							N <- 15
						} else {
							if(identical(location, "ankle")) {
								N <- 34
							} else {
								N <- 35
							}
						}
						
						for(subject in seq_len(N)) {
							system_cmd <- paste0("R --vanilla --args ", subject, " ", setting, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization, " ", update_trans,
												 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/runMethods/type/sasakiLOSOcrossval.R ")
							
							system(system_cmd, intern = TRUE)
						}
					}
				}
			}
		}
	}
}												
