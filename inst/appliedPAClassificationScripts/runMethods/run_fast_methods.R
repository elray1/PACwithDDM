#options(warn = 2, error = recover)

for(data_set in c("Mannini", "SasakiFreeLiving", "SasakiLab")) {
	for(location in c("ankle", "hip", "wrist")) {
        if(data_set %in% c("SasakiFreeLiving", "SasakiLab")) {
            class_var_levels <- c("y_category3", "y_category5", "y_intensity")
        } else {
            class_var_levels <- c("y_category4", "y_intensity")
        }
        for(class_var in class_var_levels) {
			for(fit_method in c("RF", "RFHMM", "normalHMM")) {
				setwd(file.path("C:/Stat/HMM/HMMEnsembles/HMMapplication/results", data_set, location, class_var, fit_method))

				for(reduced_trans_mat_parameterization in FALSE) {
					update_trans_levels <- TRUE
					for(update_trans in update_trans_levels) {
						
						if(identical(data_set, "SasakiFreeLiving")) {
							N <- 15
						} else if(identical(data_set, "SasakiLab")) {
							if(identical(location, "ankle")) {
								N <- 34
							} else {
								N <- 35
							}
						} else if(identical(data_set, "Mannini")) {
                            N <- 33
                        }
						
						for(subject in seq_len(N)) {
							system_cmd <- paste0("R --vanilla --args ", subject, " ", setting, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization, " ", update_trans,
												 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/runMethods/LOSO_crossval.R ")
							
							system(system_cmd, intern = TRUE)
						}
					}
				}
			}
		}
	}
}
