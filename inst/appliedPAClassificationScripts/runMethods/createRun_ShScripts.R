for(data_set in c("Mannini", "SasakiLab", "SasakiFreeLiving")) {
	if(identical(data_set, "Mannini")) {
		location_levels <- c("ankle", "wrist")
		class_var_levels <- c("y_category4", "y_intensity")
	} else {
		location_levels <- c("ankle", "wrist")
		class_var_levels <- c("y_category3", "y_intensity")
	}
	
	for(location in location_levels) {
		for(class_var in class_var_levels) {
			for(fit_method in c("parametricBoostCRF", "parametricBoostMLR")) {
				cores_req <- "50"
				span_one_host <- TRUE
				mem_req <- "2000"
				time_req <- "270:00"
				queue_req <- "long"
				
				files_path <- file.path("/home", "er71a", "HMMapplication", data_set, location, class_var, fit_method)
				setwd(files_path)
				
				basefilename <- "run_method_ShScript_reducedtime"
				lsfoutfilename <- "lsfOutput_reducedtime_"
				Routbasefilename <- "LOSOcrossval_reducedtime"
				basefilename <- paste0(basefilename, "FullTrans")
				lsfoutfilename <- paste0(lsfoutfilename, "FullTrans")
				Routbasefilename <- paste0(Routbasefilename, "FullTrans")
				
				basefilename <- paste0(basefilename, "_subject")
				lsfoutfilename <- paste0(lsfoutfilename, ".out")
				Routbasefilename <- paste0(Routbasefilename, "_subject")
				
				requestCmds <- "#!/bin/bash\n"
				
				requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n")
				if(span_one_host) {
					requestCmds <- paste0(requestCmds, "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n")
				}
				
				requestCmds <- paste0(requestCmds, "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
					"#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
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
					filename <- paste0(files_path, "/", basefilename, subject, ".sh")
					
					cat(requestCmds, file = filename)
					cat("module load R/3.2.2\n", file = filename, append = TRUE)
					cat(paste0("R CMD BATCH --vanilla \'--args ", subject, " ", data_set, " ", location, " ",
						class_var, " ", fit_method,
						"\' /home/er71a/HMMapplication/LOSO_crossval.R ", Routbasefilename, subject, ".Rout"),
						file = filename, append = TRUE)
					
					bsubCmd <- paste0("bsub < ", filename)
					
					system(bsubCmd)
				}
			}
		}
	}
}
