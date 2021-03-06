for(setting in c("freeliving")) {
	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
#		for(class_var in c("y_category3", "y_category5")) {
		for(class_var in c("y_category3")) {
#			for(fit_method in c("RFCRF", "parametricBoostCRF", "baggedFeatureSubsetGradientTreeBoostCRF")) {
#			for(fit_method in c("RFCRF")) {
			for(fit_method in c("gradientTreeBoostCRF")) {
				# set memory and processig time requirements for each method
				if(identical(fit_method, "RFCRF")) {
					cores_req <- "10"
					span_one_host <- TRUE
					mem_req <- "12000"
					time_req <- "4:00"
					queue_req <- "short"
				} else if(identical(fit_method, "RFCRFseqbag")) {
					cores_req <- "10"
					span_one_host <- TRUE
					mem_req <- "12000"
					time_req <- "4:00"
					queue_req <- "short"
				} else if(identical(fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
					cores_req <- "10"
					span_one_host <- TRUE
					mem_req <- "20000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "parametricBoostCRF")) {
					cores_req <- "20"
					span_one_host <- FALSE
					mem_req <- "1000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "L2RegularizedCRF")) {
					cores_req <- "10"
					span_one_host <- FALSE
					mem_req <- "5000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "gradientTreeBoostCRF")) {
					cores_req <- "10"
					mem_req <- "15000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "normalHMM")) {
					cores_req <- "1"
					span_one_host <- FALSE
					mem_req <- "8000"
					time_req <- "4:00"
					queue_req <- "short"
				}
				
				for(reduced_trans_mat_parameterization in c(FALSE)) {
					if(reduced_trans_mat_parameterization) {
						update_trans_levels <- TRUE
					} else {
						update_trans_levels <- TRUE
					}
					for(update_trans in update_trans_levels) {
						files_path <- file.path("/home", "er71a", "HMMapplication", "Sasaki", "intensity", setting, "2stage", "stage1", location, class_var, fit_method)
						setwd(files_path)
						
						basefilename <- "run_method_ShScript"
						lsfoutfilename <- "lsfOutput_"
						Routbasefilename <- "sasakiLOSOcrossval"
						if(reduced_trans_mat_parameterization) {
							basefilename <- paste0(basefilename, "ReducedTrans")
							lsfoutfilename <- paste0(lsfoutfilename, "ReducedTrans")
							Routbasefilename <- paste0(Routbasefilename, "ReducedTrans")
						} else {
							basefilename <- paste0(basefilename, "FullTrans")
							lsfoutfilename <- paste0(lsfoutfilename, "FullTrans")
							Routbasefilename <- paste0(Routbasefilename, "FullTrans")
						}
						if(!update_trans) {
							basefilename <- paste0(basefilename, "NoUpdate")
							lsfoutfilename <- paste0(lsfoutfilename, "NoUpdate")
							Routbasefilename <- paste0(Routbasefilename, "NoUpdate")
						}
						basefilename <- paste0(basefilename, "_subject")
						lsfoutfilename <- paste0(lsfoutfilename, ".out")
						Routbasefilename <- paste0(Routbasefilename, "_subject")
						
						requestCmds <- "#!/bin/bash\n"

						requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n")
						requestCmds <- paste0(requestCmds, "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n")
						requestCmds <- paste0(requestCmds, "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n")

						requestCmds <- paste0(requestCmds, "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
							"#BSUB -W ", time_req, " # run time\n",
							"#BSUB -q ", queue_req, " # which queue we want to run in\n")
						
						if(identical(setting, "freeliving")) {
							N <- 15
						}
						
						for(subject1 in seq_len(N - 1)) {
							for(subject2 in seq(from = subject1 + 1, to = N)) {
#						for(subject1 in c(1, 2)) {
#							for(subject2 in seq(from = subject1 + 1, to = 3)) {
								filename <- paste0(files_path, "/", basefilename, subject1, "_", subject2, ".sh")
							
								cat(requestCmds, file = filename)
								cat("module load R/3.0.2\n", file = filename, append = TRUE)
								cat(paste0("R CMD BATCH --vanilla \'--args ", subject1, " ", subject2, " ", setting, " ", location, " ", class_var, " ", fit_method, " ", reduced_trans_mat_parameterization, " ", update_trans,
													 "\' /home/er71a/HMMapplication/Sasaki/intensity/sasakiLOSOcrossval_2stage_stage1.R ", Routbasefilename, subject1, "_", subject2, ".Rout"),
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
}												
