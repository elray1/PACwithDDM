rm(list = ls())
library("rstream")
library("PACwithDDM")

# In this file, we create files with the random number generation stream starting point
# for each combination of data set, accelerometer location, classification variable, and
# fit method
save_location <- file.path(find.package("PACwithDDM"), "simStudyScripts", "rngstreams")

set.seed(81861) # this was randomly generated

# a new rstream object
rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

num_sims_per_scenario <- 100L

# represent data generation plus all possible levels for method
all_fit_methods <- c("data_gen", "parametricBoostCRF", "parametricBoostMLR", "RF", "normalHMM")

# number of substreams used by each method per simulation index
substreams_used <- list(data_gen = 1,
	parametricBoostCRF = 1001,
	parametricBoostMLR = 1001,
	RF = 1,
	normalHMM = 1)

for(obs_dist in c("obsDistNormal", "obsDistNonNormal")) {
  for(time_dependence in c("timeDependence", "noTimeDependence")) {
    for(fit_method in all_fit_methods) {
      for(sim_ind in seq_len(num_sims_per_scenario)) {
        rstream.packed(rngstream) <- TRUE
        
        stream_filename <- paste("rngstream_simStudy", obs_dist, time_dependence,
          "sim", sim_ind, sep = "_")
        save(rngstream, file = paste0(save_location, "/", fit_method, "/", stream_filename, ".rdata"))
        
        rstream.packed(rngstream) <- FALSE
        
        for(i in seq_len(eval(substreams_used[[fit_method]]))) {
          rstream.nextsubstream(rngstream)
        }
      } # sim_ind
    } # fit_method
  } # time_dependence
} # obs_dist_complex
