extract_features_one_window <- function(window_data, sampling_freq = 80) {
  # window_data is a w * sampling_freq by 3 numeric matrix:
  #   w is the window length in seconds and sampling_freq is the frequency in Hz.
  #   each row contains observations for one time point
  #   each column contains observations for one axis
  # sampling_freq is the frequency at which data were collected -- e.g., sampling_freq = 80 if data were collected at 80 Hz.
  
  # construct vector magnitude and angle
  window_data <- cbind(window_data[, c("x", "y", "z")], calc_spherical_coordinates(window_data[, c("x", "y", "z")]))
  colnames(window_data) <- c("x", "y", "z", "vm", "theta", "phi")
  
  feature_functions <- list(col_means,
    col_quantiles,
    col_lag_1_acs,
    vm_entropy,
    col_fft_features
  )
  
  feature_function_args <- list(list(window_data = quote(window_data), cols = 1:6), # col_means
    list(window_data = quote(window_data), cols = 1:6, quantile_ps = c(0.1, 0.25, 0.5, 0.75, 0.9)), # col_quantiles
    list(window_data = quote(window_data), cols = 1:4), # col_lag_1_acs
    list(window_data = quote(window_data)), # vm_entropy
    list(window_data = quote(window_data), cols = 1:4, sampling_freq = sampling_freq) # col_fft_features
  )
  
  result <- unlist(lapply(seq_along(feature_functions), function(function_ind) {
    do.call(feature_functions[[function_ind]], feature_function_args[[function_ind]])
  }))
  
  if(any(is.na(result)))
  	browser()
  
  return(result)
}

calc_spherical_coordinates <- function(window_data) {
  vm <- apply(window_data, 1, function(d1t) sqrt(sum(d1t^2)))
  theta <- acos(window_data[, 3] / vm)
  phi <- atan2(window_data[, 2], window_data[, 1])
  
  return(cbind(vm, theta, phi))
}

col_means <- function(window_data, cols) {
  temp <- colMeans(window_data[, cols, drop = FALSE])
  names(temp) <- paste0(colnames(window_data)[cols], "_mean")
  return(temp)
}

col_quantiles <- function(window_data, cols, quantile_ps) {
  temp <- as.vector(apply(window_data[, cols, drop = FALSE], 2, quantile, probs = quantile_ps))
  names(temp) <- paste0(rep(colnames(window_data[cols]), each = length(temp) / length(cols)),
    "_quantile_",
    rep(as.character(quantile_ps), times = length(cols)))
  return(temp)
}

col_lag_1_acs <- function(window_data, cols) {
  temp <- apply(window_data[, cols, drop = FALSE], 2, function(data_one_col) {
    acf(data_one_col, lag.max = 1, type = "covariance", plot = FALSE)$acf[2, 1, 1]
  })
  names(temp) <- paste0(colnames(window_data)[cols], "_lag_1_ac")
  return(temp)
}

vm_entropy <- function(window_data) {
  phat <- hist(window_data[, 4], plot = FALSE, breaks = 10)$density
  phat <- phat / sum(phat)
  phat <- phat[phat > 0]
  temp <- c("vm_entropy" = -sum(phat * log(phat)))
  return(temp)
}

col_fft_features <- function(window_data, cols, sampling_freq) {
  temp <- unlist(lapply(cols, function(col) calc_fft_features_one_col(window_data[, col], sampling_freq = sampling_freq)))
  names(temp) <- paste(rep(colnames(window_data[cols]), each = length(temp) / length(cols)),
    names(temp),
    sep = "_")
  return(temp)
}

calc_fft_features_one_col <- function(signal, sampling_freq) {
  # signal is a vector of observed values for one axis (or vector magnitude) in one window
  # sampling_freq is the frequency at which data were collected -- e.g., sampling_freq = 80 if data were collected at 80 Hz.
  
  mods <- Mod(fft(signal))[-1]
  
  n <- floor(length(mods) / 2)
  all_inds <- seq_len(n)
  
  freq <- sampling_freq * (all_inds) / (2 * n)
  mods <- mods[all_inds]
  
  total_power <- sqrt(sum(mods^2) / n)
  
  if(isTRUE(all.equal(mods, rep(0, length(mods))))) {
		like_ps <- rep(1 / length(mods), length(mods))
	} else {
		like_ps <- mods / sum(mods)
  }
  temp <- like_ps * log(like_ps)
  temp[like_ps == 0] <- 0
  entropy_spectral_density <- -sum(temp)

  dom_ind <- which.max(mods)
  dom_freq <- freq[dom_ind]
  pow_at_dom_freq <- mods[dom_ind] / sqrt(n)
  
  junk_mod <- mods[-dom_ind]
  junk_freq <- freq[-dom_ind]
  
  sec_dom_freq <- junk_freq[which.max(junk_mod)]
  pow_at_sec_dom_freq <- max(junk_mod) / sqrt(n)
  
#  like_cdf <- c(0, cumsum(mods) / sum(mods))
#  n <- length(like_cdf)
#  
#  like_quantiles <- approx(like_cdf, c(0, freq), xout = ps)$y
  
  limited_bw_inds <- all_inds[(freq >= 0.3) & (freq < 3)]
  temp <- freq[limited_bw_inds]
  limited_bw_dom_freq <- temp[which.max(mods[limited_bw_inds])]
  pow_at_limited_bw_dom_freq <- mods[freq == limited_bw_dom_freq] / sqrt(n)

	if(pow_at_dom_freq > 0) {
		relative_pow_limited_bw_dom_freq_to_overall_dom_freq <- pow_at_limited_bw_dom_freq / pow_at_dom_freq
	} else {
		relative_pow_limited_bw_dom_freq_to_overall_dom_freq <- 1
	}

  
  features <- c(total_power = total_power,
    entropy_spectral_density = entropy_spectral_density,
    dom_freq = dom_freq,
    pow_at_dom_freq = pow_at_dom_freq,
    sec_dom_freq = sec_dom_freq,
    pow_at_sec_dom_freq = pow_at_sec_dom_freq,
    limited_bw_dom_freq = limited_bw_dom_freq,
    pow_at_limited_bw_dom_freq = pow_at_limited_bw_dom_freq,
    relative_pow_limited_bw_dom_freq_to_overall_dom_freq = relative_pow_limited_bw_dom_freq_to_overall_dom_freq)
  
  return(features)	
}
