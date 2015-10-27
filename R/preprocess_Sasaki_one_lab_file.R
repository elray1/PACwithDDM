Sasaki_preprocess_one_lab_file <- function(subj, actigraph_file_path, DO_file_path, DO_adjustment_file_path, oxy_file_path, demographics_file_path, sampling_freq, window_length, drop_unlabeled = TRUE) {
  unprocessed_data <- Sasaki_combine_lab_actigraph_and_DO_files(subj, actigraph_file_path, DO_file_path, DO_adjustment_file_path, oxy_file_path, demographics_file_path, drop_unlabeled = drop_unlabeled)
  
  num_windows <- ceiling(nrow(unprocessed_data) / (window_length * sampling_freq))
  window_inds <- rep(seq_len(num_windows), each = window_length * sampling_freq)[seq_len(nrow(unprocessed_data))]
  obs_each_window <- table(window_inds)
  if(obs_each_window[length(obs_each_window)] < window_length * sampling_freq * 0.5) {
	window_inds[window_inds == window_inds[length(window_inds)]] <- window_inds[length(window_inds)] - 1
  }

  X <- tapply(seq_len(nrow(unprocessed_data)),
  									 window_inds,
  									 function(inds) {
  									 	extract_features_one_window(window_data = unprocessed_data[inds, ], sampling_freq = 80)
  									 })
  colnames_X <- names(X[[1]])
  X <- rbind.fill.matrix(lapply(X, function(comp) matrix(comp, nrow = 1)))
  colnames(X) <- colnames_X
  
  
  ## behavior
  
  y_behavior <- unname(tapply(as.character(unprocessed_data$behavior),
  	window_inds,
  	function(vals) {
  		uv <- unique(vals)
  		if(length(uv) == 1L) {
  			return(uv)
  		} else {
  			return("transition")
  		}
  	}))
  
  # if 3 or more adjacent windows have behavior type "transition", replace the behavior type of the center window(s) with "MovInter"
  temp <- rle(as.vector(y_behavior == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  rle_inds_to_modify <- which(temp$lengths >= 3 & temp$values)
  orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
  	seq(from = orig_run_start_inds[rleind] + 1, length = temp$lengths[rleind] - 2, by = 1)
  }))
  y_behavior[orig_inds_to_modify] <- "MovInter"
  
  # if 1 or 2 adjacent windows of behavior type "transition" have windows of type "MovInter" on both sides of them, set their types to "MovInter"
  # handle the sequence start and end by converting to MovInter if the condition is met on the inside
  temp <- rle(as.vector(y_behavior == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  rle_inds_to_modify <- which(temp$values &
  	c("MovInter", y_behavior)[orig_run_start_inds] == "MovInter" &
  	c(y_behavior, "MovInter")[orig_run_start_inds + temp$lengths] == "MovInter"
  )
  orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
  	seq(from = orig_run_start_inds[rleind], length = temp$lengths[rleind], by = 1)
  }))
  y_behavior[orig_inds_to_modify] <- "MovInter"
  
  # if 2 adjacent windows of behavior type "transition" are preceeded or followed by a window of type "MovInter",
  # set the type of the "transition" window next to the "MovInter" window to "MovInter"
  temp <- rle(as.vector(y_behavior == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  
  # modify left
  rle_inds_to_modify <- which(temp$values &
  	temp$length == 2 &
		c("*", y_behavior)[orig_run_start_inds] == "MovInter"
  )
  orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify]
  y_behavior[orig_inds_to_modify] <- "MovInter"
  
  # modify right
  rle_inds_to_modify <- which(temp$values &
		temp$length == 2 &
		c(y_behavior, "*")[orig_run_start_inds + 2] == "MovInter"
  )
  orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify] + 1
  y_behavior[orig_inds_to_modify] <- "MovInter"
  
  # convert to factor
  y_behavior <- factor(y_behavior, levels = c(levels(unprocessed_data$behavior), "transition"))
  
  
  ## category5
  y_category5 <- unname(tapply(as.character(unprocessed_data$category5),
  														window_inds,
  														function(vals) {
  															uv <- unique(vals)
  															if(length(uv) == 1L) {
  																return(uv)
  															} else {
  																return("transition")
  															}
  														}))
  
  # if 3 or more adjacent windows have category5 type "transition", replace the category5 type of the center window(s) with "MovInter"
  temp <- rle(as.vector(y_category5 == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  rle_inds_to_modify <- which(temp$lengths >= 3 & temp$values)
  orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
  	seq(from = orig_run_start_inds[rleind] + 1, length = temp$lengths[rleind] - 2, by = 1)
  }))
  y_category5[orig_inds_to_modify] <- "MovInter"
  
  # if 1 or 2 adjacent windows of category5 type "transition" have windows of type "MovInter" on both sides of them, set their types to "MovInter"
  # handle the sequence start and end by converting to MovInter if the condition is met on the inside
  temp <- rle(as.vector(y_category5 == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  rle_inds_to_modify <- which(temp$values &
  															c("MovInter", y_category5)[orig_run_start_inds] == "MovInter" &
  															c(y_category5, "MovInter")[orig_run_start_inds + temp$lengths] == "MovInter"
  )
  orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
  	seq(from = orig_run_start_inds[rleind], length = temp$lengths[rleind], by = 1)
  }))
  y_category5[orig_inds_to_modify] <- "MovInter"
  
  # if 2 adjacent windows of category5 type "transition" are preceeded or followed by a window of type "MovInter",
  # set the type of the "transition" window next to the "MovInter" window to "MovInter"
  temp <- rle(as.vector(y_category5 == "transition"))
  orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
  
  # modify left
  rle_inds_to_modify <- which(temp$values &
  															temp$length == 2 &
  															c("*", y_category5)[orig_run_start_inds] == "MovInter"
  )
  orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify]
  y_category5[orig_inds_to_modify] <- "MovInter"
  
  # modify right
  rle_inds_to_modify <- which(temp$values &
  															temp$length == 2 &
  															c(y_category5, "*")[orig_run_start_inds + 2] == "MovInter"
  )
  orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify] + 1
  y_category5[orig_inds_to_modify] <- "MovInter"
  
  # convert to factor
  y_category5 <- factor(y_category5, levels = c(levels(unprocessed_data$category5), "transition"))
  
  
   ## category3
   
   y_category3 <- unname(tapply(as.character(unprocessed_data$category3),
   														window_inds,
   														function(vals) {
   															uv <- unique(vals)
   															if(length(uv) == 1L) {
   																return(uv)
   															} else {
   																return("transition")
   															}
   														}))
   
   # if 3 or more adjacent windows have category3 type "transition", replace the category3 type of the center window(s) with "MovInter"
   temp <- rle(as.vector(y_category3 == "transition"))
   orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
   rle_inds_to_modify <- which(temp$lengths >= 3 & temp$values)
   orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
   	seq(from = orig_run_start_inds[rleind] + 1, length = temp$lengths[rleind] - 2, by = 1)
   }))
   y_category3[orig_inds_to_modify] <- "MovInter"
   
   # if 1 or 2 adjacent windows of category3 type "transition" have windows of type "MovInter" on both sides of them, set their types to "MovInter"
   # handle the sequence start and end by converting to MovInter if the condition is met on the inside
   temp <- rle(as.vector(y_category3 == "transition"))
   orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
   rle_inds_to_modify <- which(temp$values &
   															c("MovInter", y_category3)[orig_run_start_inds] == "MovInter" &
   															c(y_category3, "MovInter")[orig_run_start_inds + temp$lengths] == "MovInter"
   )
   orig_inds_to_modify <- unlist(lapply(rle_inds_to_modify, function(rleind) {
   	seq(from = orig_run_start_inds[rleind], length = temp$lengths[rleind], by = 1)
   }))
   y_category3[orig_inds_to_modify] <- "MovInter"
   
   # if 2 adjacent windows of category3 type "transition" are preceeded or followed by a window of type "MovInter",
   # set the type of the "transition" window next to the "MovInter" window to "MovInter"
   temp <- rle(as.vector(y_category3 == "transition"))
   orig_run_start_inds <- cumsum(c(1, temp$lengths))[seq_along(temp$lengths)]
   
   # modify left
   rle_inds_to_modify <- which(temp$values &
   															temp$length == 2 &
   															c("*", y_category3)[orig_run_start_inds] == "MovInter"
   )
   orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify]
   y_category3[orig_inds_to_modify] <- "MovInter"
   
   # modify right
   rle_inds_to_modify <- which(temp$values &
   															temp$length == 2 &
   															c(y_category3, "*")[orig_run_start_inds + 2] == "MovInter"
   )
   orig_inds_to_modify <- orig_run_start_inds[rle_inds_to_modify] + 1
   y_category3[orig_inds_to_modify] <- "MovInter"
   
   # convert to factor
   y_category3 <- factor(y_category3, levels = c(levels(unprocessed_data$category3), "transition"))
  
  
  ## intensity
  
  y_intensity <- unname(tapply(as.character(unprocessed_data$intensity),
  														window_inds,
  														function(vals) {
  															uv <- unique(vals)
  															if(length(uv) == 1L) {
  																return(uv)
  															} else {
  																return("transition")
  															}
  														}))
    
  # convert to factor
  y_intensity <- factor(y_intensity, levels = c(levels(unprocessed_data$intensity), "transition"))
  
  
  results <- list(X = X,
  	y_behavior = y_behavior,
  	y_category5 = y_category5,
  	y_category3 = y_category3,
  	y_intensity = y_intensity)
  
  return(results)
}

Sasaki_combine_lab_actigraph_and_DO_files <- function(subj, actigraph_file_path, DO_file_path, DO_adjustment_file_path, oxy_file_path, demographics_file_path, drop_unlabeled = TRUE) {
  # Read in the DO data
  temp_DO_data <- read.csv(file = file(DO_file_path), stringsAsFactors = FALSE)
  temp_DO_data$subject <- substring(temp_DO_data[, 1], 4, 5)
  temp_DO_data <- temp_DO_data[temp_DO_data$subject == subj, ]
  if(identical(subj, "23")) {
	temp_DO_data$Start[1:3] <- c("9:26:00", "9:26:31", "9:27:01")
	temp_DO_data$Stop[1:3] <- c("9:26:30", "9:27:00", "9:27:30")
  }
  
  times <- c(as.POSIXlt(strptime(temp_DO_data$Start, "%H:%M:%S")), as.POSIXlt(strptime(temp_DO_data$Stop, "%H:%M:%S")))[order(rep(seq_len(nrow(temp_DO_data)), times = 2))]
  raw_DO_data <- data.frame(Rel.Time = as.numeric(times - times[1]),
  	Subject = subj,
	Behavior = as.vector(t(cbind(temp_DO_data$Activity, rep("Unlabeled", nrow(temp_DO_data))))))

  # remove unlabeled indices of length 1 sec
  to_remove <- which(times[seq_len(length(times) - 1)] + 1 == times[1 + seq_len(length(times) - 1)])
  if(length(to_remove) > 0) {
	  times <- times[-to_remove]
	  raw_DO_data <- raw_DO_data[-to_remove, ]
  }

  length_DO <- dim(raw_DO_data)[1]
  
  DO_start_time <- times[1]
  
  raw_DO_data$time <- DO_start_time + raw_DO_data$Rel.Time

  # Add in the intensity level information, based on steady state minutes 3 - 5; METs use denominator of 3.5 (RMR is unavailable for 3 subjects)
  demographics_data <- read.csv(demographics_file_path)
  subject_mass <- demographics_data$Weight..Kg.[demographics_data$Subject.ID == paste0("ACE", subj)]

  oxy_data <- readLines(oxy_file_path)
  oxy_data <- oxy_data[sapply(oxy_data, function(comp) {substr(comp, 1, 18) != "General info field" && substr(comp, 1, 14) != "Begin artefact" && substr(comp, 1, 9) != "Occlusion" && substr(comp, 1, 13) != "Lactate entry"})]

  oxy_reltimes <- sapply(seq(from = 11, to = length(oxy_data)), function(ind) {
	temp <- as.numeric(strsplit(substr(oxy_data[ind], 1, 7), ":")[[1]])
	if(length(temp) == 2) {
		return(temp[1] * 60 + temp[2])
	} else if(length(temp) == 3) {
		return(temp[1] * 60 * 60 + temp[2] * 60 + temp[3])
	} else {
		stop("Unhandled oxy reltime case")
	}
  })
  oxy_vo <- sapply(seq(from = 11, to = length(oxy_data)), function(ind) {
	return(as.numeric(substr(oxy_data[ind], 8, 14)))
  })

  raw_DO_data$measured_METs <- NA
  raw_DO_data$intensity <- "Unlabeled"
  labeled_inds <- which(raw_DO_data$Behavior != "Unlabeled")
  raw_DO_data$intensity[labeled_inds] <- sapply(labeled_inds, function(ind) {
	if(raw_DO_data$Rel.Time[ind + 1] - raw_DO_data$Rel.Time[ind] < 5 * 60) { # not enough data to estimate METs from steady state VO2
		if(raw_DO_data$Behavior[ind] %in% c("Lying down", "Seated", "Crossword puzzles", "Playing cards", "Standing")) {
			return("Sedentary")
		} else if(raw_DO_data$Behavior[ind] %in% c("Dusting", "Gardening", "Self-care (miscellaneous)", "Laundry", "Tai-chi")) {
			return("Light")
		} else if(raw_DO_data$Behavior[ind] %in% c("Slow walk (1.8 mph)", "400m walk", "Carrying groceries", "Vacuuming", "Organizing the room", "Simulated Bowling")) {
			return("Moderate")
		} else {
			stop("Unhandled behavior type case")
		}
	} else { # 5 minutes or more of data -- estimate METs from steady state VO2
		start_time <- raw_DO_data$time[ind]
		end_time <- raw_DO_data$time[ind + 1]
		temp_ind <- which(as.POSIXlt(strptime(temp_DO_data$Start, "%H:%M:%S")) == start_time)

		# start 2 minutes in
		temp <- as.numeric(strsplit(temp_DO_data$Start.Oxy[temp_ind], ":")[[1]])
		oxy_start_reltime <- temp[1] * 60 * 60 + temp[2] * 60 + temp[3] + 2 * 60

		oxy_inds_keep <- (oxy_reltimes >= oxy_start_reltime & oxy_reltimes < oxy_start_reltime + 3 * 60)
		raw_DO_data$measured_METs[ind] <<- mean(oxy_vo[oxy_inds_keep], na.rm = TRUE) / subject_mass / 3.5
		if(raw_DO_data$measured_METs[ind] < 1.5) {
			return("Sedentary")
		} else if(raw_DO_data$measured_METs[ind] >= 1.5 && raw_DO_data$measured_METs[ind] < 3.0) {
			return("Light")
		} else if(raw_DO_data$measured_METs[ind] >= 3.0 && raw_DO_data$measured_METs[ind] <= 6.0) {
			return("Moderate")
		} else {
			return("Vigorous")
		}
	}
  })

  # Read in the DO adjustment data and make adjustments to times, if specified
  if(!missing(DO_adjustment_file_path)) {
    # perform Ray adjustments
    temp <- read.csv(file = file(DO_adjustment_file_path), skip = 0, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
    overall_adjustment <- as.numeric(temp[1, 2])
    
    temp <- read.csv(file = file(DO_adjustment_file_path), skip = 1, header = TRUE, stringsAsFactors = FALSE)
    temp$Time <- as.POSIXlt(strptime(temp[, 1], "%H:%M:%S"))
    if(!isTRUE(all.equal(as.numeric(temp$Time), as.numeric(times))))
      stop("Error: Require that relative times in DO file match original relative times in DO adjustment file.")
    
    raw_DO_data$Rel.Time <- raw_DO_data$Rel.Time + overall_adjustment + as.numeric(temp[, 2])
    
    raw_DO_data$time <- DO_start_time + raw_DO_data$Rel.Time
  }


  # Read in the accelerometer data
  raw_accel_data <- read.csv(actigraph_file_path, as.is = TRUE, skip = 10, header = FALSE, stringsAsFactors = FALSE)
  
  temp <- read.csv(actigraph_file_path, skip = 0, nrows = 10, header = FALSE)
  accel_start_time_row <- which(apply(temp, 1, function(line) grepl("Start Time*", line, perl = TRUE)))
  accel_start_time <- substring(temp[accel_start_time_row, 1], 12)
  accel_start_time <- as.POSIXlt(strptime(accel_start_time,"%H:%M:%S"))
  raw_accel_data$time <- accel_start_time + seq(from = 0, length = dim(raw_accel_data)[1], by = 1/80)
  
  
  # Combine the DO data and the accelerometer data.
  # This has two steps:
  # (1) Reduce the data range of the raw.accel.data to only that which is covered by DO observations
  # (2) Add the DO activity description at each time point to the accelerometer data
  
  # Step (1): reduce data range.
  
  # first_ind is the first index of the raw.accel.data matrix which falls within the range of times covered by direct observation.
  # specifically, it is the smallest index i such that raw_accel_data$time[i] is >= DO_start_time
  first_ind <- min(which(raw_accel_data$time >= raw_DO_data$time[1]))
  
  # last_ind is the last index of the raw_accel_data matrix which falls within the range of times covered by direct observation.
  # specifically, it is the largest index i such that raw.accel.times[i] is <= raw_DO_data$time[nrow(raw_DO_data)]
  last_ind <- max(which(raw_accel_data$time <= raw_DO_data$time[nrow(raw_DO_data)]))
  
  combined_data <- as.data.frame(raw_accel_data[seq(from = first_ind, to = last_ind), c(2, 3, 1, 4)])
  names(combined_data) <- c("x", "y", "z", "time")
  combined_data$rel_time <- seq(from = 0, length = nrow(combined_data), by = 1/80)
  
  # corr_ind is a vector of the length of the combined_data dataframe that gives the corresponding row in the DO data for each row in the combined_data
  corr_ind <- apply(matrix(combined_data$time), 1, function(tt) {
    temp <- which(raw_DO_data$time[seq(from = 1, to = nrow(raw_DO_data) - 1)] <= tt &
      raw_DO_data$time[seq(from = 2, to = nrow(raw_DO_data))] >= tt)
    if(length(temp) == 0)
      stop(paste0("Error: unable to find corresponding ind for ", tt))
    return(min(temp))
  })

  # copy the DO classification of the behavior and intensity for each time point in the accelerometer data.
  combined_data$behavior <- raw_DO_data[corr_ind, "Behavior"]
  combined_data$intensity <- raw_DO_data[corr_ind, "intensity"]

  
  # A little cleanup:
  #  - Drop unlabeled time if requested
  #  - Group classifications
  #  - Turn classifications into factors
  #  - MET values
  
  if(drop_unlabeled) {
	combined_data <- combined_data[combined_data$behavior != "Unlabeled", ]
  }

  combined_data$category5 <- NA
  combined_data$category5[combined_data$behavior %in% c("Standing")] <- "Standing"
  combined_data$category5[combined_data$behavior %in% c("Lying down", "Seated", "Crossword puzzles", "Playing cards")] <- "Sedentary"
  combined_data$category5[combined_data$behavior %in% c("Dusting", "Gardening", "Vacuuming", "Self-care (miscellaneous)", "Laundry", "Organizing the room")] <- "MovInter"
  combined_data$category5[combined_data$behavior %in% c("Slow walk (1.8 mph)", "400m walk", "Carrying groceries")] <- "Locomotion"
  combined_data$category5[combined_data$behavior %in% c("Tai-chi", "Simulated Bowling")] <- "Recreational"
#  if(any(!(combined_data$category5 %in% c("Standing", "Sedentary", "MovInter", "Locomotion", "Recreational"))))
#    stop("Error in setting category5")

  combined_data$category3 <- as.character(combined_data$category5)
  combined_data$category3[combined_data$category3 %in% c("Sedentary", "Standing")] <- "Sedentary/Standing"
  combined_data$category3[combined_data$category3 %in% c("MovInter", "Recreational")] <- "MovInter"


	combined_data$behavior <- factor(combined_data$behavior, levels = sort(unique(combined_data$behavior)))
	combined_data$category5 <- factor(combined_data$category5, levels = c("Standing", "Sedentary", "MovInter", "Locomotion", "Recreational"))
	combined_data$category3 <- factor(combined_data$category3, levels = c("Sedentary/Standing", "MovInter", "Locomotion"))
	if(drop_unlabeled) {
		combined_data$intensity <- factor(combined_data$intensity, levels = c("Sedentary", "Light", "Moderate"))
	} else {
		combined_data$intensity <- factor(combined_data$intensity, levels = c("Sedentary", "Light", "Moderate", "Unlabeled"))
	}

  return(combined_data)
}
