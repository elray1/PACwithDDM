rm(list = ls())

library("R.matlab")
library("rayThesis")

rayThesis_location <- find.package("rayThesis")

temp <- readMat(file.path(rayThesis_location, "extdata", "Mannini", "SD2010_wrist_12800_correctionON_90Hz.mat"))

subject <- temp$Labels["Participant", , 1][[1]]

# use just the features they computed already
wrist_features <- data.frame(temp$Data["P", ,][[1]])

posture_labels <- unlist(temp$Info[, , 1][, 1]$posture.labels)
y_posture <- factor(posture_labels[as.integer(as.vector(temp$Labels["Lab1w", , 1][[1]]))])

activity_labels <- unlist(temp$Info[, , 1][, 1]$activity.labels)
y_behavior <- factor(activity_labels[as.integer(as.vector(temp$Labels["Lab2w", , 1][[1]]))])

category4_labels <- c("ambulation", "cycling", "other", "sedentary")
y_category4 <- factor(category4_labels[as.integer(as.vector(temp$Labels["Lab5w", , 1][[1]]) - 1)], levels = category4_labels)

y_intensity <- rep(NA, length(y_posture))
y_intensity[y_behavior %in% c("none", "unknown ", "sitting-at-computer:-internet-search ", "sitting-at-computer:-typing ", "sitting:-reading ", "sitting:-writing ")] <- "Sedentary"
y_intensity[y_behavior %in% c("sorting-files-or-paperwork")] <- "Light"
y_intensity[y_behavior %in% c("sweeping-with-broom", "painting:-roller ", "painting:-brush ", "stairs:-inside-and-down ",
	"lifting-10-pound-box ", "cycling:-70-rpm_-50-watts_-.7-kg ", "cycling:-outdoor-downhill ", "cycling:-outdoor-level ",
	"treadmill:-3.0-mph_-0%", "treadmill:-2.0-mph_-0%", "treadmill:-2.0-mph_-9%", "treadmill:-4.0-mph_-0%",
	"walking:-natural", "walking:-prescribed", "carrying-load ")] <- "Moderate"
y_intensity[y_behavior %in% c("cycling:-outdoor-uphill ", "stairs:-inside-and-up", "treadmill:-3.0-mph_-6%", "treadmill:-3.0-mph_-9%", "treadmill:-6.0-mph_-0%")] <- "Vigorous"

y_intensity <- factor(y_intensity, levels = c("Sedentary", "Light", "Moderate", "Vigorous"))

ManniniWristCorrected <- lapply(unique(subject), function(subj) {
	inds <- which(subject == subj)

	return(list(X = as.matrix(wrist_features[inds, ]),
		y_behavior = y_behavior[inds],
		y_posture = y_posture[inds],
		y_category4 = y_category4[inds],
		y_intensity = y_intensity[inds],
		subject = subj,
		location = "wrist"))
})

save(ManniniWristCorrected, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/ManniniWristCorrected.rdata")




temp <- readMat(file.path(rayThesis_location, "extdata", "Mannini", "SD2010_ankle_12800_correctionON_90Hz.mat"))

subject <- temp$Labels["Participant", , 1][[1]]

# use just the features they computed already
wrist_features <- data.frame(temp$Data["P", ,][[1]])

posture_labels <- unlist(temp$Info[, , 1][, 1]$posture.labels)
y_posture <- factor(posture_labels[as.integer(as.vector(temp$Labels["Lab1w", , 1][[1]]))])

activity_labels <- unlist(temp$Info[, , 1][, 1]$activity.labels)
y_behavior <- factor(activity_labels[as.integer(as.vector(temp$Labels["Lab2w", , 1][[1]]))])

category4_labels <- c("ambulation", "cycling", "other", "sedentary")
y_category4 <- factor(category4_labels[as.integer(as.vector(temp$Labels["Lab5w", , 1][[1]]) - 1)], levels = category4_labels)

y_intensity <- rep(NA, length(y_posture))
y_intensity[y_behavior %in% c("none", "unknown ", "sitting-at-computer:-internet-search ", "sitting-at-computer:-typing ", "sitting:-reading ", "sitting:-writing ")] <- "Sedentary"
y_intensity[y_behavior %in% c("sorting-files-or-paperwork")] <- "Light"
y_intensity[y_behavior %in% c("sweeping-with-broom", "painting:-roller ", "painting:-brush ", "stairs:-inside-and-down ",
	"lifting-10-pound-box ", "cycling:-70-rpm_-50-watts_-.7-kg ", "cycling:-outdoor-downhill ", "cycling:-outdoor-level ",
	"treadmill:-3.0-mph_-0%", "treadmill:-2.0-mph_-0%", "treadmill:-2.0-mph_-9%", "treadmill:-4.0-mph_-0%",
	"walking:-natural", "walking:-prescribed", "carrying-load ")] <- "Moderate"
y_intensity[y_behavior %in% c("cycling:-outdoor-uphill ", "stairs:-inside-and-up", "treadmill:-3.0-mph_-6%", "treadmill:-3.0-mph_-9%", "treadmill:-6.0-mph_-0%")] <- "Vigorous"

y_intensity <- factor(y_intensity, levels = c("Sedentary", "Light", "Moderate", "Vigorous"))

ManniniAnkleCorrected <- lapply(unique(subject), function(subj) {
	inds <- which(subject == subj)

	return(list(X = as.matrix(wrist_features[inds, ]),
		y_behavior = y_behavior[inds],
		y_posture = y_posture[inds],
		y_category4 = y_category4[inds],
		y_intensity = y_intensity[inds],
		subject = subj,
		location = "wrist"))
})

save(ManniniAnkleCorrected, file = "C:/Stat/HMM/HMMEnsembles/rayThesis/data/ManniniAnkleCorrected.rdata")

