rm(list = ls())
library("rstream")

set.seed(416455) # this was randomly generated

# a new rstream object
rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))

# dimension of observations for Sasaki data
D <- 77

# Build a data frame with the number of substreams used for each combination,
# including methods that are not run by this code, but that do use rng substreams.
# The substream offset is then 1 + the total number of substreams used by all earlier methods

# possible levels for method
all_fit_methods <- c("RFCRF", "RFCRFseqbag", "parametricBoostCRF", "RF", "normalHMM",
	"gradientTreeBoostCRF", "baggedFeatureSubsetGradientTreeBoostCRF", "MLRHMM", "RFHMM")

# create the data frame discussed above -- initialize with all combinations, then remove irrelevant combinations
substream_df <- data.frame(fit_method = rep(all_fit_methods, each = 4),
	case = rep(c("1stage", "1stageusingtrueclass", "2stage_stage1", "2stage_stage2"), times = length(all_fit_methods)),
	stringsAsFactors = FALSE)

# set number of substreams used per subject -- depends on fit_method, case, D, and number of classes
substream_df$substreams_used <- 0

substream_df$substreams_used[substream_df$fit_method == "RFCRF"] <- 1001
substream_df$substreams_used[substream_df$fit_method == "RFCRFseqbag"] <- 1001
substream_df$substreams_used[substream_df$fit_method == "parametricBoostCRF"] <- 1001
substream_df$substreams_used[substream_df$fit_method == "RF"] <- 1
substream_df$substreams_used[substream_df$fit_method == "normalHMM"] <- 1
substream_df$substreams_used[substream_df$fit_method == "gradientTreeBoostCRF" & substream_df$case %in% c("1stage", "2stage_stage1")] <- (2 + 10 * D * 30 + 1000)
substream_df$substreams_used[substream_df$fit_method == "gradientTreeBoostCRF" & substream_df$case %in% c("1stageusingtrueclass", "2stage_stage2")] <- (2 + 10 * (D + 4) * 30 + 1000)
substream_df$substreams_used[substream_df$fit_method == "baggedFeatureSubsetGradientTreeBoostCRF" & substream_df$case %in% c("1stage", "2stage_stage1")] <- (2 + 10 * D * 30 + 1000)
substream_df$substreams_used[substream_df$fit_method == "baggedFeatureSubsetGradientTreeBoostCRF" & substream_df$case %in% c("1stageusingtrueclass", "2stage_stage2")] <- (2 + 10 * (D + 4) * 30 + 1000)

substream_df <- rbind(c("", "", 0), substream_df)


for(setting in c("freeliving", "lab")) {
#for(setting in c("freeliving")) {
	# number of subjects
	if(identical(setting, "freeliving")) {
		N <- 15
	} else {
		N <- 35
	}

	for(location in c("ankle", "hip", "wrist")) {
#	for(location in c("ankle")) {
		# create rstream object with a new seed.  For programming simplicity and reduced computation time,
		# we use a new seed for each combination of setting, location, and class_var,
		# with the same seed value (but different substreams used) for levels of
		# fit_method, reduced_trans_mat_parameterization, and update_trans within
		# the combination of setting, location, and class_var.
		# The seed values are randomly generated.

		# Get rng substream offset corresponding to the combination of fit_method, reduced_trans_mat_parameterization, and update_trans:

		for(current_scenario_row_ind in seq(from = 2, to = nrow(substream_df))) {
			if(current_scenario_row_ind == 2 && !(identical(setting, "freeliving") && identical(location, "ankle"))) {
				# prevent overlap with baggedFeatureSubsetGradientTreeBoostCRF from previous scenario
				for(i in seq_len(2 + 10 * (D + 4) * 30 + 1000)) {
					rstream.nextsubstream(rngstream)
				}
			}

			if(identical(substream_df$case[current_scenario_row_ind], "2stage_stage1")) {
				for(subj in seq_len(N - 1)) {
					for(subj2 in seq(from = subj + 1, to = N)) {
						if(subj == 1 && subj2 == 2) {
							# prevent overlap with last subject, previous fit method
							for(i in seq_len(substream_df$substreams_used[current_scenario_row_ind - 1])) {
								rstream.nextsubstream(rngstream)
							}
						} else {
							for(i in seq_len(substream_df$substreams_used[current_scenario_row_ind])) {
								rstream.nextsubstream(rngstream)
							}
						}
						rstream.packed(rngstream) <- TRUE
						stream_filename <- paste("rngstream_Sasaki", setting, location, substream_df$fit_method[current_scenario_row_ind],
							"case", substream_df$case[current_scenario_row_ind], "subj", subj, subj2, sep = "_")
						save(rngstream, file = paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/intensity/rngstreams/", stream_filename, ".rdata"))
						rstream.packed(rngstream) <- FALSE
					} # subj2
				} # subj
			} else {
				for(subj in seq_len(N)) {
					if(subj == 1) {
						# prevent overlap with last subject, previous fit method
						for(i in seq_len(substream_df$substreams_used[current_scenario_row_ind - 1])) {
							rstream.nextsubstream(rngstream)
						}
					} else {
						for(i in seq_len(substream_df$substreams_used[current_scenario_row_ind])) {
							rstream.nextsubstream(rngstream)
						}
					}
					rstream.packed(rngstream) <- TRUE
					stream_filename <- paste("rngstream_Sasaki", setting, location, substream_df$fit_method[current_scenario_row_ind],
						"case", substream_df$case[current_scenario_row_ind], "subj", subj, sep = "_")
					save(rngstream, file = paste0("C:/Stat/HMM/HMMEnsembles/rayThesis/inst/appliedPAClassificationScripts/Sasaki/intensity/rngstreams/", stream_filename, ".rdata"))
					rstream.packed(rngstream) <- FALSE
				} # subj
			} # case "2stage_stage1"
		} # current_scenario_row_ind
	} # location
} # setting
