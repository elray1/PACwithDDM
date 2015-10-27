for(setting in c("lab")) {
	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
		for(truth in c(FALSE, TRUE)) {
			for(fit_method in c("RFCRF", "RFCRFseqbag", "baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "gradientTreeBoostCRF")) {
#			for(fit_method in c("RFCRFseqbag", "baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "gradientTreeBoostCRF")) {
#			for(fit_method in c("RFCRF", "RFCRFseqbag")) {
				# set memory and processing time requirements for each method
				if(identical(fit_method, "RFCRF")) {
					cores_req <- "10"
					mem_req <- "12000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "RFCRFseqbag")) {
					cores_req <- "10"
					mem_req <- "12000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
					cores_req <- "10"
					mem_req <- "20000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "parametricBoostCRF")) {
					cores_req <- "10"
					mem_req <- "8000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "L2RegularizedCRF")) {
					cores_req <- "10"
					mem_req <- "5000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "gradientTreeBoostCRF")) {
					cores_req <- "10"
					mem_req <- "15000"
					time_req <- "48:00"
					queue_req <- "long"
				}

				if(truth) {
					truth_folder <- "truth"
				} else {
					truth_folder <- "notruth"
				}
				
				files_path <- file.path("/home", "er71a", "HMMapplication", "Sasaki", "intensity", setting, "1stage", truth_folder, location, fit_method)
				setwd(files_path)
						
				basefilename <- "run_method_ShScript_subject"
				lsfoutfilename <- "lsfOutput.out"
				Routbasefilename <- "sasakiLOSOcrossval_subject"
				
				requestCmds <- "#!/bin/bash\n"

				requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n")
				requestCmds <- paste0(requestCmds, "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n")
				requestCmds <- paste0(requestCmds, "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n")
				requestCmds <- paste0(requestCmds, "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
					"#BSUB -W ", time_req, " # run time\n",
					"#BSUB -q ", queue_req, " # which queue we want to run in\n")
						
				if(identical(setting, "freeliving")) {
					N <- 15
				} else if(identical(setting, "lab")) {
					if(identical(location, "ankle")) {
						N <- 34
					} else {
						N <- 35
					}
				}
				
#				for(subject in seq_len(N)) {
				for(subject in c(1, seq(from = 3, to = N))) {
#				for(subject in 2) {
					filename <- paste0(files_path, "/", basefilename, subject, ".sh")
							
					cat(requestCmds, file = filename)
					cat("module load R/3.0.2\n", file = filename, append = TRUE)
					cat(paste0("R CMD BATCH --vanilla \'--args ", subject, " ", setting, " ", location, " ", fit_method, " ", truth,
										 "\' /home/er71a/HMMapplication/Sasaki/intensity/sasakiLOSOcrossval_1stage.R ", Routbasefilename, subject, ".Rout"),
							file = filename, append = TRUE)
							
					bsubCmd <- paste0("bsub < ", filename)
							
					system(bsubCmd)
				}
			}
		}
	}
}												
