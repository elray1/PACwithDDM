rm(list = ls())
library("rstream")

set.seed(81861) # this was randomly generated

# a new rstream object
rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

num_sims_per_scenario <- 100L

# represent data generation plus all possible levels for method
all_fit_methods <- c("data_gen", "RFCRF", "RFCRFseqbag", "parametricBoostCRF", "RF", "normalHMM", "unregularizedParametricCRF", "L2RegularizedCRF",
	"gradientTreeBoostCRF", "baggedFeatureSubsetGradientTreeBoostCRF", "MLRHMM", "RFHMM")

# number of substreams used by each method per simulation index
substreams_used <- list(data_gen = 1,
	RFCRF = 1001,
	RFCRFseqbag = 1001,
	parametricBoostCRF = 1001,
	RF = 1,
	normalHMM = 1,
	unregularizedParametricCRF = 1,
	L2RegularizedCRF = 11,
	gradientTreeBoostCRF = expression(2 + 10 * D * 30 + 1000),
	baggedFeatureSubsetGradientTreeBoostCRF = expression(2 + 10 * D * 30 + 1000),
	MLRHMM = 1,
	RFHMM = 0)

D <- 50

for(obs_dist_complex in c(FALSE, TRUE)) {
for(bayes_error_rate_high in c(FALSE, TRUE)) {
for(redundant_features_informative in c(FALSE, TRUE)) {
for(fit_method in all_fit_methods) {
	if(obs_dist_complex) {
		obs_dist <- "obsDistNormal"
	} else {
		obs_dist <- "obsDistNonNormal"
	}

	if(bayes_error_rate_high) {
		bayes_error <- "BayesErrorLarge"
	} else {
		bayes_error <- "BayesErrorSmall"
	}

	if(redundant_features_informative) {
		redundant_features <- "redundantFeaturesInformative"
	} else {
		redundant_features <- "redundantFeaturesNonInformative"
	}

	for(sim_ind in seq_len(num_sims_per_scenario)) {
		rstream.packed(rngstream) <- TRUE

		stream_filename <- paste("rngstream_simStudy", obs_dist, bayes_error, redundant_features,
			"sim", sim_ind, sep = "_")
		save(rngstream, file = paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/simStudyScripts/rngstreams/", fit_method, "/", stream_filename, ".rdata"))

		rstream.packed(rngstream) <- FALSE

		for(i in seq_len(eval(substreams_used[[fit_method]]))) {
			rstream.nextsubstream(rngstream)
		}
	} # subj
} # fit_method
} # redundant_features_informative
} # bayes_error_rate_high
} # obs_dist_complex
