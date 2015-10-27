options(warn = 2)

for(setting in c("freeliving")) {
	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
#		for(class_var in c("y_category3", "y_category5")) {
		for(class_var in c("y_category3")) {
			for(fit_method in c("RF", "RFHMM")) {
#			for(fit_method in c("RFHMM")) {
				setwd(file.path("C:/Stat/HMM/HMMEnsembles/HMMapplication/Sasaki/intensity/results", setting, "2stage", "stage1", location, class_var, fit_method))

				if(fit_method %in% c("RFHMM", "MLRHMM", "normalHMM")) {
#					reduced_trans_mat_parameterization_levels <- c(TRUE, FALSE)
#					reduced_trans_mat_parameterization_levels <- TRUE
					reduced_trans_mat_parameterization_levels <- FALSE
				} else {
					reduced_trans_mat_parameterization_levels <- TRUE
				}
				for(reduced_trans_mat_parameterization in reduced_trans_mat_parameterization_levels) {
					update_trans_levels <- TRUE
					for(update_trans in update_trans_levels) {
						
						if(identical(setting, "freeliving")) {
							N <- 15
						}
						
						for(subject1 in seq(from = 1, to = N - 1)) {
							for(subject2 in seq(from = subject1 + 1, to = N)) {
#						for(subject1 in c(1, 2)) {
#							for(subject2 in seq(from = subject1 + 1, to = 3)) {
								system_cmd <- paste0("R --vanilla --args ", subject1, " ", subject2, " ", setting, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization, " ", update_trans,
												 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/runMethods/intensity/sasakiLOSOcrossval_2stage_stage1.R ")
							
								system(system_cmd, intern = TRUE)
							}
						}
					}
				}
			}
		}
	}
}
