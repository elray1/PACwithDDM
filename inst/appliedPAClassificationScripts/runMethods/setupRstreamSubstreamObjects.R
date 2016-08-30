rm(list = ls())
library("rstream")
library("PACwithDDM")

# In this file, we create files with the random number generation stream starting point
# for each combination of data set, accelerometer location, classification variable, and
# fit method
save_location <- file.path(find.package("PACwithDDM"), "appliedPAClassificationScripts", "rngstreams")

# First, build a data frame with the number of substreams used for each fit_method

# create the data frame discussed above -- initialize with all combinations, then remove irrelevant combinations
substream_df <- data.frame(fit_method = c("RF", "normalHMM", "parametricBoostCRF", "parametricBoostMLR", "RFHMM"),
	substreams_used = c(1, # RF: 1 stream for all subjects
	1, # normalHMM: 1 stream for all subjects
	10001, # parametricBoostCRF: 1 to generate sequence bag groups + M_bag = 10000, one for each bag group
	10001, # parametricBoostMLR: 1 to generate sequence bag groups + M_bag = 10000, one for each bag group
	0), # normalFMM: the method does not require random number generation
	stringsAsFactors = FALSE)



# create rstream object.  the seed values are randomly generated.
set.seed(327073)
rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

# Iterate over combinations of data set, accelerometer location, classification variable, fit method,
# transition matrix parameterization, and subject.  For each combination, save the rstream object and then advance
# by the number of streams used by the given fit_method.
for(data_set in c("Mannini", "SasakiFreeLiving", "SasakiLab")) {
	# set values for number of subjects, locations, and classification variables,
    # depending on data_set
    if(identical(data_set, "Mannini")) {
        N <- 33
        locations <- c("ankle", "wrist")
        class_vars <- c("y_category4", "y_intensity")
    } else if(identical(data_set, "SasakiFreeLiving")) {
        N <- 15
        locations <- c("ankle", "wrist")
        class_vars <- c("y_category3", "y_intensity")
    } else {
    		N <- 35
        locations <- c("ankle", "wrist")
        class_vars <- c("y_category3", "y_intensity")
    }

	for(location in locations) {
		for(class_var in class_vars) {
			for(current_scenario_row_ind in seq_len(nrow(substream_df))) {
				for(subj in seq_len(N)) {
					rstream.packed(rngstream) <- TRUE
					stream_filename <- paste0(paste("rngstream", data_set, location, class_var,
                        substream_df$fit_method[current_scenario_row_ind],
												"subj", subj, sep = "_"), ".rdata")
					save(rngstream, file = file.path(save_location, stream_filename))
					rstream.packed(rngstream) <- FALSE
                    for(i in seq_len(substream_df$substreams_used[current_scenario_row_ind])) {
                        rstream.nextsubstream(rngstream)
                    }
                } # subj
			} # current_scenario_row_ind
		} # class_var
	} # location
} # data_set
