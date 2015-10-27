rm(list = ls())
library("rstream")
library("DDMforPA")

# In this file, we create files with the random number generation stream starting point
# for each combination of data set, accelerometer location, classification variable, and
# fit method
save_location <- file.path(find.package("DDMforPA"), "appliedPAClassificationScripts", "rngstreams")

# First, build a data frame with the number of substreams used for each fit_method and
# whether or not a reduced transition matrix is relevant for that method

# possible levels for method
all_fit_methods <- c("RF", "normalHMM", "L2RegularizedCRF", "RFHMM")

# whether the reduced_trans_mat_parameterization is relevant to each method
reduced_trans_mat_parameterization_relevant <- c(FALSE, TRUE, TRUE, TRUE)

# create the data frame discussed above -- initialize with all combinations, then remove irrelevant combinations
substream_df <- data.frame(fit_method = rep(all_fit_methods, each = 2),
	reduced_trans = rep(c(TRUE, FALSE), times = length(all_fit_methods)),
	stringsAsFactors = FALSE)

# remove rows where reduced_trans_mat_parameterization is not relevant and reduced_trans is FALSE
to_remove <- substream_df$fit_method %in% all_fit_methods[!reduced_trans_mat_parameterization_relevant] & !substream_df$reduced_trans
substream_df <- substream_df[!to_remove, ]

# set number of substreams used per subject -- depends on fit_method only
num_substreams_used <- c("RF"=1, # RF: 1 stream for all subjects
	"normalHMM"=1, # normalHMM: 1 stream for all subjects in each case of FullTrans and ReducedTrans
	"L2RegularizedCRF"=11, # L2RegularizedCRF: 1 to generate sequence crossvalidation groups + K_crossval = 10, one for each crossvalidation group
	"RFHMM"=0) # RFHMM: the method uses the fits created for the RF method above.

substream_df$substreams_used <- sapply(substream_df$fit_method, function(fm) {
	num_substreams_used[all_fit_methods == fm]
})


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
        locations <- c("ankle", "hip", "wrist")
        class_vars <- c("y_category4", "y_intensity")
    } else if(identical(data_set, "SasakiFreeLiving")) {
        N <- 15
        locations <- c("ankle", "hip", "wrist")
        class_vars <- c("y_category5", "y_category3", "y_intensity")
    } else {
		N <- 35
        locations <- c("ankle", "hip", "wrist")
        class_vars <- c("y_category5", "y_category3", "y_intensity")
    }

	for(location in locations) {
		for(class_var in class_vars) {
			for(current_scenario_row_ind in seq_len(nrow(substream_df))) {
				for(subj in seq_len(N)) {
					rstream.packed(rngstream) <- TRUE
					stream_filename <- paste0(paste("rngstream", data_set, location, class_var,
                        substream_df$fit_method[current_scenario_row_ind],
						"reducedTrans", as.character(substream_df$reduced_trans[current_scenario_row_ind]),
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
