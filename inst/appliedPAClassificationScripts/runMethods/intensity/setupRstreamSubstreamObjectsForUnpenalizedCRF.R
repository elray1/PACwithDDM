rm(list = ls())
library("rstream")

# dimension of observations for Sasaki data
D <- 77

# load the last saved rngstream offset object
prev_setting <- "lab"
prev_location <- "wrist"
prev_fit_method <- "RFHMM"
prev_case <- "2stage_stage2"
prev_subj <- 35L
prev_stream_filename <- paste("rngstream_Sasaki", prev_setting, prev_location, prev_fit_method, "case", prev_case, "subj", prev_subj, sep = "_")

load(file = paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/intensity/rngstreams/", prev_stream_filename, ".rdata"))
rstream.packed(rngstream) <- FALSE

# advance the corresponding number of substreams: 25302
for(i in seq_len(25302)) {
	rstream.nextsubstream(rngstream)
}

# possible levels for method
new_fit_method <- "unpenalizedParametricCRF"

# number of substreams used per subject: 1
num_substreams_used <- 1


for(setting in c("freeliving", "lab")) {
	# number of subjects
	if(identical(setting, "freeliving")) {
		N <- 15
	} else {
		N <- 35
	}

	for(location in c("ankle", "hip", "wrist")) {
		# create rstream object with a new seed.  For programming simplicity and reduced computation time,
		# we use a new seed for each combination of setting, location, and class_var,
		# with the same seed value (but different substreams used) for levels of
		# fit_method, reduced_trans_mat_parameterization, and update_trans within
		# the combination of setting, location, and class_var.
		# The seed values are randomly generated.

		# Get rng substream offset corresponding to the combination of fit_method, reduced_trans_mat_parameterization, and update_trans:
		for(case in c("1stage", "1stageusingtrueclass")) {
			for(subj in seq_len(N)) {
				rstream.packed(rngstream) <- TRUE
				stream_filename <- paste("rngstream_Sasaki", setting, location, new_fit_method,
					"case", case, "subj", subj, sep = "_")
				save(rngstream, file = paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/intensity/rngstreams/", stream_filename, ".rdata"))
				save(rngstream, file = paste0("C:/Program Files/R/R-3.0.2/library/rayThesis/appliedPAClassificationScripts/Sasaki/intensity/rngstreams/", stream_filename, ".rdata"))
				rstream.packed(rngstream) <- FALSE
				for(i in seq_len(num_substreams_used)) {
					rstream.nextsubstream(rngstream)
				}
			} # subj
		} # truth
	} # location
} # setting
