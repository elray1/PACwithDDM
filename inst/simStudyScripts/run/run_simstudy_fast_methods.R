rm(list = ls())
options(warn = 2, error = recover)

for(obs_dist_normal in c(FALSE, TRUE)) {
	for(bayes_error_rate_high in c(FALSE, TRUE)) {
#	for(bayes_error_rate_high in c(FALSE)) {
#	for(bayes_error_rate_high in c(TRUE)) {
		for(redundant_features_informative in c(FALSE, TRUE)) {
			for(fit_method in c("RF", "RFHMM", "MLRHMM", "normalHMM", "unregularizedParametricCRF", "BayesRuleTestParams")) {
#			for(fit_method in c("RF", "RFHMM", "MLRHMM", "unregularizedParametricCRF", "BayesRuleTestParams")) {
#			for(fit_method in c("normalHMM")) {
					sim_inds <- seq_len(50)
#					sim_inds <- 2
					
					for(sim_ind in sim_inds) {
						system_cmd <- paste0("R --vanilla --args ", sim_ind, " ", obs_dist_normal, " ", redundant_features_informative, " ", bayes_error_rate_high, " ", fit_method,
											 " < C:/Stat/HMM/HMMEnsembles/rayThesis/inst/simStudyScripts/run/simStudyOneCase.R")
						
						system(system_cmd, intern = TRUE)
					}
				}
			}
		}
	}
}
