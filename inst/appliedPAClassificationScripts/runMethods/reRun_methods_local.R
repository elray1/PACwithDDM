options(warn = 2, error = recover)

for(data_set in c("Mannini", "SasakiLab", "SasakiFreeLiving")) {
    if(identical(data_set, "Mannini")) {
		location_levels <- c("ankle", "wrist")
	} else {
		location_levels <- c("ankle", "hip", "wrist")
	}
	for(location in location_levels) {
		if(identical(data_set, "Mannini")) {
			class_var_levels <- c("y_category4", "y_intensity")
		} else {
			class_var_levels <- c("y_category3", "y_category5", "y_intensity")
		}

		for(class_var in class_var_levels) {
			for(fit_method in c("L2RegularizedCRF")) {
				for(reduced_trans_mat_parameterization in c(FALSE)) {
#					files_path <- file.path("/home", "er71a", "HMMapplication", data_set, location, class_var, fit_method)
					files_path <- file.path("C:", "Stat", "HMM", "PACwithDDM", "inst", "results", data_set, location, class_var, fit_method)
#				    files_path <- file.path("F:", "Evan", "PACwithDDM-linux", "pacwithddm", "inst", "results", data_set, location, class_var, fit_method)
				  setwd(files_path)
					
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
#					    if(!file.exists(file.path("F:", "Evan", "PACwithDDM-linux", "pacwithddm", "inst", "results", data_set, location, class_var, fit_method, paste0("results_FullTrans_subject", subject, ".Rdata")))) {
						if(!file.exists(file.path("C:", "Stat", "HMM", "PACwithDDM", "inst", "results", data_set, location, class_var, fit_method, paste0("results_FullTrans_subject", subject, ".Rdata")))) {
							system_cmd <- paste0("R --vanilla --args ", subject, " ", data_set, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization,
#							    " < F:/Evan/PACwithDDM-linux/pacwithddm/inst/appliedPAClassificationScripts/runMethods/LOSO_crossval.R")
								" < C:/Stat/HMM/PACwithDDM/inst/appliedPAClassificationScripts/runMethods/LOSO_crossval.R")

							system(system_cmd, intern = TRUE)
						} else {
							cat(paste0("skipping subject ", subject, "\n"))
						}
					}
				}
			}
		}
	}
}
