for(obs_dist_normal in c(FALSE, TRUE)) {
	for(bayes_error_rate_high in c(FALSE, TRUE)) {
	#for(bayes_error_rate_high in c(FALSE)) {
	#for(bayes_error_rate_high in c(TRUE)) {
#		for(redundant_features_informative in c(FALSE, TRUE)) {
		for(redundant_features_informative in TRUE) {
			for(fit_method in c("RFCRF", "RFCRFseqbag", "baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "gradientTreeBoostCRF")) {
#			for(fit_method in c("RFCRFseqbag", "baggedFeatureSubsetGradientTreeBoostCRF", "parametricBoostCRF", "gradientTreeBoostCRF")) {
#			for(fit_method in c("gradientTreeBoostCRF")) {
#			for(fit_method in c("RFCRF")) {
				# set memory and processig time requirements for each method
				if(identical(fit_method, "RFCRF")) {
					cores_req <- "10"
					mem_req <- "12000"
					time_req <- "4:00"
					queue_req <- "short"
				} else if(identical(fit_method, "RFCRFseqbag")) {
					cores_req <- "10"
					mem_req <- "12000"
					time_req <- "4:00"
					queue_req <- "short"
				} else if(identical(fit_method, "baggedFeatureSubsetGradientTreeBoostCRF")) {
					cores_req <- "10"
					mem_req <- "20000"
					time_req <- "48:00"
					queue_req <- "long"
				} else if(identical(fit_method, "parametricBoostCRF")) {
					cores_req <- "10"
					mem_req <- "2000"
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
			
				if(obs_dist_normal) {
				        obs_dist_folder <- "obsDistNormal"
				} else {
				        obs_dist_folder <- "obsDistNonNormal"
				}
				
				if(bayes_error_rate_high) {
				        bayes_error_folder <- "BayesErrorLarge"
				} else {
				        bayes_error_folder <- "BayesErrorSmall"
				}
				
				if(redundant_features_informative) {
				        redundant_features_folder <- "redundantFeaturesInformative"
				} else {
				        redundant_features_folder <- "redundantFeaturesNonInformative"
				}

				files_path <- file.path("/home", "er71a", "simStudy", "logs", obs_dist_folder, bayes_error_folder, redundant_features_folder, fit_method)
				setwd(files_path)
				
				basefilename <- "run_method_ShScript_simind"
				lsfoutfilename <- "lsfOutput.out"
				Routbasefilename <- "simStudy_simind"
			
				requestCmds <- "#!/bin/bash\n"

				requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n")
				requestCmds <- paste0(requestCmds, "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n")
				requestCmds <- paste0(requestCmds, "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n")
				requestCmds <- paste0(requestCmds, "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
					"#BSUB -W ", time_req, " # run time\n",
					"#BSUB -q ", queue_req, " # which queue we want to run in\n")
					
#				sim_inds <- seq_len(50)
#				sim_inds <- 2
				sim_inds <- c(1, 3:50)
				
				for(sim_ind in sim_inds) {
					filename <- paste0(files_path, "/", basefilename, sim_ind, ".sh")
					
					cat(requestCmds, file = filename)
					cat("module load R/3.0.2\n", file = filename, append = TRUE)
					cat(paste0("R CMD BATCH --vanilla \'--args ", sim_ind, " ", obs_dist_normal, " ", redundant_features_informative, " ", bayes_error_rate_high, " ", fit_method,
									 "\' /home/er71a/simStudy/simStudyOneCase.R ", Routbasefilename, sim_ind, ".Rout"),
						file = filename, append = TRUE)
			
					bsubCmd <- paste0("bsub < ", filename)
					
					system(bsubCmd)
				}
			}
		}
	}
}
