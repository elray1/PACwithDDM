Sasaki_preprocess_one_free_living_file <- function(actigraph_file_path, DO_file_path, DO_adjustment_file_path, sampling_freq, window_length, drop_private = TRUE) {
  unprocessed_data <- Sasaki_combine_free_living_actigraph_and_DO_files(actigraph_file_path, DO_file_path, DO_adjustment_file_path, drop_private = drop_private)
  
  num_windows <- ceiling(nrow(unprocessed_data) / (window_length * sampling_freq))
  window_inds <- rep(seq_len(num_windows), each = window_length * sampling_freq)[seq_len(nrow(unprocessed_data))]
  
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
  
  y_intensity <- unname(tapply(as.character(unprocessed_data$behavior_mod_2),
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
  y_intensity <- factor(y_intensity, levels = c(levels(unprocessed_data$behavior_mod_2), "transition"))
  
  
  results <- list(X = X,
  	y_behavior = y_behavior,
  	y_category5 = y_category5,
  	y_category3 = y_category3,
  	y_intensity = y_intensity)
  
  return(results)
}

Sasaki_combine_free_living_actigraph_and_DO_files <- function(actigraph_file_path, DO_file_path, DO_adjustment_file_path, drop_private = TRUE) {
  # Read in the DO data
  raw_DO_data <- read.csv(file = file(DO_file_path), skip = 4, stringsAsFactors = FALSE)
  raw_DO_data <- rbind(raw_DO_data[raw_DO_data$Event.Type == "State start", c(1, 2, 5, 6, 7, 8)], raw_DO_data[dim(raw_DO_data)[1], c(1, 2, 5, 6, 7, 8)])
  colnames(raw_DO_data) <- c("Rel.Time", "Subject", "Behavior", "Behavior.Mod1", "Behavior.Mod2", "Event.Type")
  length_DO <- dim(raw_DO_data)[1]
  
  temp <- read.csv(file = file(DO_file_path), skip = 1, nrows = 2, header = FALSE, stringsAsFactors = FALSE)
  DO_start_time <- as.POSIXlt(strptime(temp[1, 2], "%H:%M:%S")) + (as.numeric(temp[2, 2]) / 1000)
  
  Sasaki_DO_corrections_used <- FALSE
  
  # Read in the DO adjustment data and make adjustments to times, if specified
  if(!missing(DO_adjustment_file_path)) {
    if(identical(substr(DO_adjustment_file_path, nchar(DO_adjustment_file_path) - 10, nchar(DO_adjustment_file_path)), "DO.log2.csv")) {
      # perform Sasaki adjustments -- replace DO data by values specified in Sasaki's file, but keep start time
      Sasaki_DO_corrections_used <- TRUE
      
      DO_adjustment_records <- read.csv(DO_adjustment_file_path, header = TRUE, stringsAsFactors = FALSE)
      DO_adjustment_records <- DO_adjustment_records[DO_adjustment_records$Observation.Name == raw_DO_data$Subject[1], ]
      
      raw_DO_data <- rbind(raw_DO_data[1, ],
                           data.frame(Rel.Time = DO_adjustment_records$Relative.Time..seconds,
                                Subject = raw_DO_data$Subject[1],
                                Behavior = DO_adjustment_records$Behavior,
                                Behavior.Mod1 = DO_adjustment_records$Behavior.Modifier.1,
                                Behavior.Mod2 = DO_adjustment_records$Behavior.Modifier.2,
                                Event.Type = DO_adjustment_records$Event.Type))
      
      raw_DO_data$time <- c(DO_start_time, as.POSIXct(strptime(substr(DO_adjustment_records$adj_starttimes, nchar(DO_adjustment_records$adj_starttimes) - 7, nchar(DO_adjustment_records$adj_starttimes)), "%H:%M:%S")))
    } else {
      # perform Ray adjustments
      temp <- read.csv(file = file(DO_adjustment_file_path), skip = 0, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
      overall_adjustment <- as.numeric(temp[1, 2])
      
      temp <- read.csv(file = file(DO_adjustment_file_path), skip = 1, header = TRUE, stringsAsFactors = FALSE)
      if(!isTRUE(all.equal(as.numeric(temp[, 1]), as.numeric(raw_DO_data$Rel.Time))))
        stop("Error: Require that relative times in DO file match original relative times in DO adjustment file.")
      
      raw_DO_data$Rel.Time <- raw_DO_data$Rel.Time + overall_adjustment + as.numeric(temp[, 2])
      
      raw_DO_data$time <- DO_start_time + raw_DO_data$Rel.Time
    }
  } else {
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
  first_ind <- min(which(raw_accel_data$time >= DO_start_time))
  
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

  # copy the DO classification of the behavior and behavior modifier for each time point in the accelerometer data.
  combined_data[, c("behavior", "behavior_mod_1", "behavior_mod_2")] <- raw_DO_data[corr_ind, c("Behavior", "Behavior.Mod1", "Behavior.Mod2")]
  
  # A little cleanup:
  #  - Remove entries with behavior == "private" if requested
  #  - Rename Sedent to Sedentary.
  #  - Fill in missing intensity classifications based on CPA values corresponding to recorded behavior
  #  - Group classifications
  #  - Turn classifications into factors
  if(drop_private)
	  combined_data <- combined_data[combined_data$behavior != "Private", ]
  combined_data$behavior_mod_2[combined_data$behavior_mod_2 == "Sedent"] <- "Sedentary"
  combined_data$behavior_mod_2[combined_data$behavior_mod_2 == "" & combined_data$behavior %in% c("Sitting", "Standing")] <- "Sedentary" # CPA values of 1.3 METs
  combined_data$behavior_mod_2[combined_data$behavior_mod_2 == "" & combined_data$behavior %in% c("Driving")] <- "Light" # CPA value of 2.5 METs
  
  combined_data$category5 <- NA
  combined_data$category5[combined_data$behavior %in% c("Standing", "StUpper")] <- "Standing"
  combined_data$category5[combined_data$behavior %in% c("Lying", "Sitting", "SitUpper", "Driving")] <- "Sedentary"
  combined_data$category5[combined_data$behavior %in% c("MovInter")] <- "MovInter"
  combined_data$category5[combined_data$behavior %in% c("Walking", "WalkLoad", "WalkInc", "Stairs")] <- "Locomotion"
  combined_data$category5[combined_data$behavior %in% c("ResistEx", "BalancEx", "AerobEx")] <- "Recreational"
#  if(any(!(combined_data$category5 %in% c("Standing", "Sedentary", "Household", "Locomotion", "Recreational"))))
#    stop("Error in setting category5")
  
  combined_data$category3 <- as.character(combined_data$category5)
  combined_data$category3[combined_data$category3 %in% c("Sedentary", "Standing")] <- "Sedentary/Standing"
  combined_data$category3[combined_data$category3 %in% c("MovInter", "Recreational")] <- "MovInter"
#  if(any(!(combined_data$category3 %in% c("Sedentary/Standing", "MovingIntermittently", "Locomotion"))))
#    stop("Error in setting category3")

	combined_data$category5 <- factor(combined_data$category5, levels = c("Standing", "Sedentary", "MovInter", "Locomotion", "Recreational"))
	combined_data$category3 <- factor(combined_data$category3, levels = c("Sedentary/Standing", "MovInter", "Locomotion"))


  if(!Sasaki_DO_corrections_used) {
  	if(drop_private) {
  		combined_data$behavior <- factor(combined_data$behavior,
  			levels = c("StUpper", "SitUpper", "MovInter", "Walking", "Standing", "WalkLoad", "Stairs", "Sitting",
  				"Driving", "ResistEx", "BalancEx", "Lying", "WalkInc", "AerobEx"))
  	} else {
  		combined_data$behavior <- factor(combined_data$behavior,
	      levels = c("Private", "StUpper", "SitUpper", "MovInter", "Walking", "Standing", "WalkLoad", "Stairs", "Sitting",
 	       "Driving", "ResistEx", "BalancEx", "Lying", "WalkInc", "AerobEx"))
  	}
  } else {
    combined_data$behavior <- factor(combined_data$behavior,
                                     levels = sort(unique(combined_data$behavior)))
  }
  combined_data$behavior_mod_1 <- factor(combined_data$behavior_mod_1, levels = c("", "Indoor", "Outdoor"))
  combined_data$behavior_mod_2 <- factor(combined_data$behavior_mod_2, levels = c("", "Sedentary", "Light", "Moderate", "Vigorous"))
  
  return(combined_data)
}
