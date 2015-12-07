for(data_set in c("Mannini", "SasakiLab", "SasakiFreeLiving")) {
	if(identical(data_set, "Mannini")) {
		location_levels <- c("ankle", "wrist")
#		location_levels <- c("wrist")
	} else {
		location_levels <- c("ankle", "hip", "wrist")
	}
	for(location in location_levels) {
#	for(location in c("ankle")) {
		if(identical(data_set, "Mannini")) {
			class_var_levels <- c("y_category4", "y_intensity")
		} else {
			class_var_levels <- c("y_category3", "y_category5", "y_intensity")
		}

		for(class_var in class_var_levels) {
			for(fit_method in c("L2RegularizedCRF")) {
				cores_req <- "10"
				span_one_host <- TRUE
				mem_req <- "5000"
				time_req <- "12:00"
				queue_req <- "long"
				
				for(reduced_trans_mat_parameterization in c(FALSE)) {
					files_path <- file.path("/home", "er71a", "HMMapplication", data_set, location, class_var, fit_method)
					setwd(files_path)
					
					basefilename <- "run_method_ShScript"
					lsfoutfilename <- "lsfOutput_"
					Routbasefilename <- "LOSOcrossval"
					if(reduced_trans_mat_parameterization) {
						basefilename <- paste0(basefilename, "ReducedTrans")
						lsfoutfilename <- paste0(lsfoutfilename, "ReducedTrans")
						Routbasefilename <- paste0(Routbasefilename, "ReducedTrans")
					} else {
						basefilename <- paste0(basefilename, "FullTrans")
						lsfoutfilename <- paste0(lsfoutfilename, "FullTrans")
						Routbasefilename <- paste0(Routbasefilename, "FullTrans")
					}

					basefilename <- paste0(basefilename, "_subject")
					lsfoutfilename <- paste0(lsfoutfilename, ".out")
					Routbasefilename <- paste0(Routbasefilename, "_subject")
					
					requestCmds <- "#!/bin/bash\n"

					requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n")
					if(span_one_host) {
						requestCmds <- paste0(requestCmds, "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n")
					}

					if(fit_method %in% c("gradientTreeBoostCRF")) {
						requestCmds <- paste0(requestCmds, "#BSUB -R \"1*{rusage[mem=", mem_req_first, "]} + 9*{rusage[mem=", mem_req_others, "]}\" # ask for memory\n")
					} else {
						requestCmds <- paste0(requestCmds, "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n")
					}

					requestCmds <- paste0(requestCmds, "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
						"#BSUB -W ", time_req, " # run time\n",
						"#BSUB -q ", queue_req, " # which queue we want to run in\n")
					
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
						if(!file.exists(file.path("/home", "er71a", "HMMapplication", "results", data_set, location, class_var, fit_method, paste0("results_FullTrans_subject", subject, ".Rdata")))) {
							filename <- paste0(files_path, "/", basefilename, subject, ".sh")
							
							cat(requestCmds, file = filename)
							cat("module load R/3.0.2\n", file = filename, append = TRUE)
							cat(paste0("R CMD BATCH --vanilla \'--args ", subject, " ", data_set, " ", location, " ",
								class_var, " ", fit_method, " ", reduced_trans_mat_parameterization,
								"\' /home/er71a/HMMapplication/LOSO_crossval.R ", Routbasefilename, subject, ".Rout"),
								file = filename, append = TRUE)
							
							bsubCmd <- paste0("bsub < ", filename)
							
							system(bsubCmd)
						}
					}
				}
			}
		}
	}
}												
