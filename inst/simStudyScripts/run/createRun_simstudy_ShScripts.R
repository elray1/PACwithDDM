for(obs_dist_normal in c(FALSE, TRUE)) {
  for(time_dependence in c(FALSE, TRUE)) {
    for(fit_method in c("parametricBoostCRF", "parametricBoostMLR")) {
      # set memory and processing time requirements for each method
      cores_req <- "16"
      mem_req <- "6000"
      time_req <- "360:00"
      queue_req <- "long"
      
      if(obs_dist_normal) {
        obs_dist_folder <- "obsDistNormal"
      } else {
        obs_dist_folder <- "obsDistNonNormal"
      }
      
      if(time_dependence) {
        time_dependence_folder <- "timeDependence"
      } else {
        time_dependence_folder <- "noTimeDependence"
      }
      
      files_path <- file.path("/home", "er71a", "simStudy", "logs", obs_dist_folder, time_dependence_folder, fit_method)
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
      
      sim_inds <- seq_len(50)
      for(sim_ind in sim_inds) {
        filename <- paste0(files_path, "/", basefilename, sim_ind, ".sh")
        
        cat(requestCmds, file = filename)
        cat("module load R/3.2.2\n", file = filename, append = TRUE)
        cat(paste0("R CMD BATCH --vanilla \'--args ", sim_ind, " ", obs_dist_normal, " ", time_dependence, " ", fit_method,
          "\' /home/er71a/simStudy/simStudyOneCase.R ", Routbasefilename, sim_ind, ".Rout"),
          file = filename, append = TRUE)
        
        bsubCmd <- paste0("bsub < ", filename)
        
        system(bsubCmd)
      }
    }
  }
}
