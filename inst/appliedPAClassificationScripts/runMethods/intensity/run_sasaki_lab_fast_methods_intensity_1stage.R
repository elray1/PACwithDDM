options(warn = 2)

for(setting in c("lab")) {
#	for(truth in c(FALSE, TRUE)) {
	for(truth in TRUE) {
		if(truth) {
			truth_folder <- "truth"
		} else {
			truth_folder <- "notruth"
		}

		for(location in c("ankle", "hip", "wrist")) {
#		for(location in c("ankle")) {
#			for(fit_method in c("RF", "RFHMM", "MLRHMM", "normalHMM", "unregularizedParametricCRF")) {
#			for(fit_method in c("RFHMM", "MLRHMM")) {
			for(fit_method in c("unregularizedParametricCRF")) {
				setwd(file.path("C:", "Stat", "HMM", "HMMEnsembles", "HMMapplication", "Sasaki", "intensity", "results", setting, "1stage", truth_folder, location, fit_method))

				if(fit_method %in% c("RFHMM", "MLRHMM", "normalHMM")) {
#					reduced_trans_mat_parameterization_levels <- c(TRUE, FALSE)
					reduced_trans_mat_parameterization_levels <- FALSE
#					reduced_trans_mat_parameterization_levels <- TRUE
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
						
						for(subject in seq_len(N)) {
#						for(subject in 2) {
							system_cmd <- paste0("R --vanilla --args ", subject, " ", setting, " ", location, " ", fit_method, " ", truth,
												 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/runMethods/intensity/sasakiLOSOcrossval_1stage.R ")
							
							system(system_cmd, intern = TRUE)
						}
					}
				}
			}
		}
	}
}
