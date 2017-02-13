rm(list = ls())
options(warn = 2, error = recover)

for(obs_dist_normal in c(FALSE, TRUE)) {
	for(time_dependence in c(FALSE, TRUE)) {
		for(fit_method in c("RF", "normalFMM", "normalHMM")) {
			sim_inds <- seq_len(50)
			
			for(sim_ind in sim_inds) {
				system_cmd <- paste0("R --vanilla --args ", sim_ind, " ", obs_dist_normal, " ", time_dependence, " ", fit_method,
					 " < C:/Stat/HMM/PACwithDDM/inst/simStudyScripts/run/simStudyOneCase.R")
				
				system(system_cmd, intern = TRUE)
			}
		}
	}
}
