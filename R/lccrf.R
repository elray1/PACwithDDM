initialize_new_CRF <- function(D = 1,
    S = 1,
    log_pi = log(rep(1/S, S)),
    log_trans_mat = matrix(log(rep(1/S, S^2)), nrow = S)) {
    if(D < 1)
        stop("D must be >= 1.")
    
    crf <- .Call("initialize_new_CRF_C", as.integer(D), as.integer(S),
        as.numeric(log_pi), as.numeric(log_trans_mat))
    
    return(crf)
}

set_CRF_pi_from_log <- function(HMM_Xptr, log_pi) {
    return(.Call("set_CRF_pi_from_log_C", HMM_Xptr, log_pi))
}

set_CRF_trans_matrix_from_log <- function(HMM_Xptr, log_trans_mat) {
    return(.Call("set_CRF_trans_matrix_from_log_C", HMM_Xptr, log_trans_mat))
}


calc_omegas_from_omega_hat <- function(omega_hat, S, crf_control) {
  if(crf_control$reduced_trans_mat_parameterization) {
	omegas <- lapply(seq_len(S - 1), function(s) {
		temp <- rep(omega_hat, S)
		temp[s] <- 0
		return(temp)
	})
	omegas[[S]] <- rep(omega_hat, S - 1)
  } else {
	  omegas <- lapply(seq_len(S - 1), function(s) {
		omega_hat[(s - 1) * (S) + seq_len(S)]
	  })
	  omegas[[S]] <- omega_hat[(S - 1) * (S) + seq_len(S - 1)]
  }

  return(omegas)
}

calc_log_trans_mat_from_omegas_CRF <- function(omegas) {
  S <- length(omegas)
  log_trans_mat <- matrix(NA, nrow = S, ncol = S)
  
  for(s in seq_len(S - 1))
    log_trans_mat[s, ] <- omegas[[s]]
  log_trans_mat[S, ] <- c(omegas[[S]], 0)

#  rowsums <- logspace_sum_matrix_rows(log_trans_mat)
#  log_trans_mat <- sweep(log_trans_mat, 1, rowsums, `-`)
  
  return(log_trans_mat)
}

calc_omegas_from_log_trans_mat_CRF <- function(log_trans_mat) {
  S <- nrow(log_trans_mat)
  
  result <- lapply(seq_len(S - 1), function(s) {
    omegas <- log_trans_mat[s, ]
  })
  result[[S]] <- log_trans_mat[s, seq_len(S - 1)]
  
  return(result)
}

set_trans_mat_from_omegas_CRF <- function(CRF_Xptr, omegas) {
  set_CRF_trans_matrix_from_log(CRF_Xptr, calc_log_trans_mat_from_omegas_CRF(omegas))
}


lccrf_compute_initial_betas_MLR <- function(CRF_Xptr, active_cols, data_concat, data_split, T_sub, S,
		omegas, prev_log_obs_probs_by_state, crf_control) {
	X <- data_concat$X[, active_cols, drop = FALSE]
	if(!is.factor(data_concat$y)) {
		y_vec <- factor(data_concat$y)
	} else {
		y_vec <- data_concat$y
	}

	set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)

	num_active_cols <- length(active_cols)

	sub_start_inds <- c(1, 1 + cumsum(T_sub)[seq_len(length(T_sub) - 1)])
  
    obs_weights <- rep(1, sum(T_sub))

	# get offsets
	offsets <- calc_marginal_class_probs_given_log_obs_probs_CRF(CRF_Xptr, data_split, prev_log_obs_probs_by_state,
		log = TRUE)

	n_obs_classes <- length(unique(data_concat$y))
	if(identical(1L, n_obs_classes)) {
		result <- rep(-10, S * num_active_cols)
		class_ind <- (levels(data_concat$y) == unique(data_concat$y))
		result[class_ind - 1 + seq_len(num_active_cols)] <- 0
	} else {
		reduced_levels <- levels(y_vec)[levels(y_vec) %in% data_concat$y]
		reduced_level_inds <- which(levels(y_vec) %in% reduced_levels)
		y_vec <- factor(y_vec, levels = reduced_levels)
		y_for_mlr <- matrix(0, length(y_vec), n_obs_classes)
		y_for_mlr[cbind(seq_along(y_vec), sapply(y_vec, function(yy) which(reduced_levels == yy)))] <- 1

		if(isTRUE(all.equal(X[, 1], rep(1, nrow(X)))) && ncol(X) > 1)
			X <- X[, 2:ncol(X), drop = FALSE]
		
		temp <- multinom(y_for_mlr ~ X + offset(offsets[, reduced_level_inds]), weights = obs_weights)

		result <- rep(-10, S * num_active_cols)

		inds <- c()
		for(rli in reduced_level_inds)
			inds <- c(inds, (rli - 1) * num_active_cols + seq_len(num_active_cols))
		
		result[inds] <- c(rep(0, num_active_cols), as.numeric(t(coef(temp))))
	}

	result <- result - rep(result[(S - 1) * num_active_cols + seq_len(num_active_cols)], S)
	result <- result[seq_len((S - 1) * num_active_cols)]

	return(result)
}





lccrf_control <- function(crf_control, fit_method = "rf", reduced_trans_mat_parameterization = TRUE,
    M_boost = 100, M_boost_search_threshold = 20,
    num_active_vars = 1, active_var_search_method = "random", add_intercept = TRUE, quadratic = FALSE, center_scale = TRUE, max_attempts = 100L, beta_initialization = "random", beta_penalty_factor = "crossval-select",
    max_tree_depth = 2L,
	bag_method = "none", M_bag = 1, sequence_bag_sample_size, sequence_bag_sample_replace = TRUE, timepoint_bag_sample_size, timepoint_bag_sample_size_proportion,
    max_rf_refit_iter, rf_refit_tol, update_log_rhos = TRUE, update_transition_matrix = TRUE,
    optim_method = "L-BFGS-B", crossval_method = "fit-all", K_crossval = 10, crossval_step_only = FALSE, parallel_method = "none",
    rng_method = "default", rng_seed = NULL, rstream_substream_offset = 0,
    save_freq = 20, save_path, save_filename_base) {
	if(missing(crf_control))
		crf_control <- list()
	if(!is.list(crf_control))
		stop("If specified, crf_control must be a list.")

	if(!missing(fit_method))
	  crf_control$fit_method <- fit_method
	if(!missing(reduced_trans_mat_parameterization))
	  crf_control$reduced_trans_mat_parameterization <- reduced_trans_mat_parameterization

	if(!missing(M_boost))
		crf_control$M_boost <- M_boost
	if(!missing(M_boost_search_threshold))
		crf_control$M_boost_search_threshold <- M_boost_search_threshold

	if(!missing(num_active_vars))
		crf_control$num_active_vars <- num_active_vars
	if(!missing(active_var_search_method))
		crf_control$active_var_search_method <- active_var_search_method

	crf_control$add_intercept <- add_intercept
	if(!missing(quadratic))
		crf_control$quadratic <- quadratic
	crf_control$center_scale <- center_scale

	crf_control$max_attempts <- max_attempts

	if(!missing(beta_initialization))
	  crf_control$beta_initialization <- beta_initialization
	crf_control$beta_penalty_factor <- beta_penalty_factor
	
	if(!missing(max_tree_depth))
	  crf_control$max_tree_depth <- max_tree_depth
  
	if(!missing(bag_method))
	  crf_control$bag_method <- bag_method
	if(!missing(M_bag))
	  crf_control$M_bag <- M_bag
	if(!missing(sequence_bag_sample_size))
	  crf_control$sequence_bag_sample_size <- sequence_bag_sample_size
	if(!missing(sequence_bag_sample_replace))
	  crf_control$sequence_bag_sample_replace <- sequence_bag_sample_replace
	if(!missing(timepoint_bag_sample_size))
	  crf_control$timepoint_bag_sample_size <- timepoint_bag_sample_size
	if(!missing(timepoint_bag_sample_size_proportion))
	  crf_control$timepoint_bag_sample_size_proportion <- timepoint_bag_sample_size_proportion
	
	if(!missing(max_rf_refit_iter))
	  crf_control$max_rf_refit_iter <- max_rf_refit_iter
	if(!missing(rf_refit_tol))
	  crf_control$rf_refit_tol <- rf_refit_tol
	if(!missing(update_log_rhos))
	  crf_control$update_log_rhos <- update_log_rhos
  crf_control$update_transition_matrix <- update_transition_matrix

	if(!missing(optim_method))
		crf_control$optim_method <- optim_method

	if(!missing(crossval_method))
		crf_control$crossval_method <- crossval_method
	if(!missing(K_crossval))
		crf_control$K_crossval <- K_crossval
	crf_control$crossval_step_only <- crossval_step_only

	if(!missing(parallel_method))
	  crf_control$parallel_method <- parallel_method
  
	if(!missing(rng_method))
	  crf_control$rng_method <- rng_method
	if(!missing(rng_seed))
	  crf_control$rng_seed <- rng_seed
    if(!missing(rstream_substream_offset))
      crf_control$rstream_substream_offset <- rstream_substream_offset
	
	if(!missing(save_freq))
		crf_control$save_freq <- save_freq
	if(!missing(save_path))
		crf_control$save_path <- save_path
	if(!missing(save_filename_base))
		crf_control$save_filename_base <- save_filename_base

	return(crf_control)
}

lccrf_validate_crf_control <- function(crf_control) {
}


lccrf <- function(data_concat, crf_control, rngstream, resume_from) {
	if(missing(crf_control)) {
		crf_control <- lccrf_crf_control()
	} else {
		lccrf_validate_crf_control(crf_control)
	}

	if(missing(resume_from))
		resume_from <- NULL
  
  # set up RNG
  if(!is.null(crf_control$rng_seed))
    set.seed(crf_control$rng_seed)
  
  if(identical(crf_control$rng_method, "rstream")) {
	if(missing(rngstream)) {
	    rngstream <- new("rstream.mrg32k3a", seed = sample(1:10000, 6, rep = FALSE))
    } else if(rstream.packed(rngstream)) {
		rstream.packed(rngstream) <- FALSE
	}
    
    for(i in seq_len(crf_control$rstream_substream_offset))
      rstream.nextsubstream(rngstream)
    
    rstream.RNG(rngstream)
  } else {
    rngstream <- NULL
  }
  
  if(identical(crf_control$parallel_method, "none")) {
    parallel_lapply_method <- lapply
  } else if(identical(crf_control$parallel_method, "snowfall")) {
    parallel_lapply_method <- sfLapply
  } else {
    stop("Invalid parallel_method: must be none or snowfall")
  }

	# sample size
	N <- length(unique(data_concat$subject))
	T_total <- nrow(data_concat$X)


	# clean up X
	if(crf_control$fit_method %in%  c("parametric-boost", "parametric-L2-penalized-MLE")) {
  		# add leading column of 1's if necessary
  		if(!isTRUE(all.equal(data_concat$X[, 1], rep(1, T_total))) && crf_control$add_intercept)
  			data_concat$X <- cbind(rep(1, T_total), data_concat$X)
  
  		# center and scale covariates if requested
		if(crf_control$center_scale) {
  			cent <- apply(data_concat$X, 2, mean)
  			cent[1] <- 0
  
  			scal <- apply(data_concat$X, 2, sd)
  			scal[scal == 0] <- 1
  
  			data_concat$X <- scale(data_concat$X, center = cent, scale = scal)
		} else {
			cent <- rep(0, ncol(data_concat$X))
			scal <- rep(1, ncol(data_concat$X))
		}

		# add trailing column of 0's for L2 penalized method -- otherwise, coefficients would not be treated the same.
		if(identical(crf_control$fit_method, "parametric-L2-penalized-MLE")) {
			data_concat$X <- cbind(data_concat$X, rep(0, nrow(data_concat$X)))
		}
	}

	# number of variables
	D <- ncol(data_concat$X)

	# clean up y
	if(!is.factor(data_concat$y)) {
		warning("y should be a factor: converting now")
		data_concat$y <- factor(data_concat$y)
	}
	orig_y_levels <- levels(data_concat$y)
	levels(data_concat$y) <- seq_along(levels(data_concat$y))

	# generate groups used for training component model parameters and selecting M
	if("sequence" %in% crf_control$bag_method) {
		# generate bags to use in bagging, if applicable
		sequence_bag_groups <- lapply(seq_len(crf_control$M_bag), function(ind) sample.int(N, size = crf_control$sequence_bag_sample_size, replace = crf_control$sequence_bag_sample_replace))
		inds_where_sampled_all <- sapply(sequence_bag_groups, function(bag_group) {
			return(length(unique(bag_group)) == N)
		})
		while(any(inds_where_sampled_all) && !(crf_control$sequence_bag_sample_size == N && !crf_control$sequence_bag_sample_replace)) {
			sequence_bag_groups[inds_where_sampled_all] <- lapply(seq_len(sum(inds_where_sampled_all)), function(ind) sample.int(N, size = crf_control$sequence_bag_sample_size, replace = crf_control$sequence_bag_sample_replace))
			inds_where_sampled_all <- sapply(sequence_bag_groups, function(bag_group) {
				return(length(unique(bag_group)) == N)
			})
		}
	} else if(identical(crf_control$bag_method, "none")) {
		# generate crossvalidation groups to use for selecting M, if applicable
		crossval_group_inds <- seq_len(crf_control$K_crossval)
		crossval_groups <- sample(rep(crossval_group_inds, times = ceiling(N / crf_control$K_crossval))[seq_len(N)])
		crossval_groups <- lapply(crossval_group_inds, function(group) which(crossval_groups == group))
	} else {
		sequence_bag_groups <- lapply(seq_len(crf_control$M_bag), function(ind) seq_len(N))
	}

  if(identical(crf_control$fit_method, "rf")) {
  	if(identical(crf_control$rng_method, "rstream"))
  	  rstream.packed(rngstream) <- TRUE
      
    temp_result <- rf_crf(data_concat = data_concat, sequence_bag_groups = sequence_bag_groups, rngstream = rngstream, rstream_offset = 0, crf_control = crf_control)

  	if(identical(crf_control$rng_method, "rstream"))
  	  rstream.packed(rngstream) <- FALSE
	# rng is not used again, or I would need to update rngstream

  } else if(crf_control$fit_method %in% c("parametric-boost", "gradient-tree-boost")) {
    # if necessary, perform cross-validation to select tree depths and/or num active vars
    # this call affects rngstream state -- note that no rng is done below without first resetting rngstream state
    if(identical(crf_control$fit_method, "gradient-tree-boost") &&
  			(identical(crf_control$max_tree_depth, "crossval-select") || identical(crf_control$num_active_vars, "crossval-select"))) {
	  	if(identical(crf_control$rng_method, "rstream"))
	 	  rstream.packed(rngstream) <- TRUE
      
  		temp <- lccrf_crossval_select_tree_depth_and_num_active_vars_gradient_tree(data_concat, crf_control, rngstream, rstream_offset = 0)
		if(crf_control$crossval_step_only) {
			return(temp)
		}

	  	if(identical(crf_control$rng_method, "rstream")) {
	 	  rstream.packed(rngstream) <- FALSE
		  for(i in seq_len(1 + crf_control$K_crossval * D * 30))
			rstream.nextsubstream(rngstream)
		}

  		if(identical(crf_control$max_tree_depth, "crossval-select")) {
  			if(identical(crf_control$num_active_vars, "crossval-select")) { # crossval-select both depth and num_active_vars
  				crf_control$max_tree_depth <- list(temp[1], "crossval-select")
  				crf_control$num_active_vars <- list(temp[2], "crossval-select")
  			} else { # crossval-select depth only
  				crf_control$max_tree_depth <- list(temp[1], "crossval-select")
  			}
  		} else { # crossval-select num_active_vars only
  			crf_control$num_active_vars <- list(temp[1], "crossval-select")
  		}
 	}
  	
  	# get final estimates
  	if(identical(crf_control$bag_method, "sequence")) {
  		if(identical(crf_control$rng_method, "rstream"))
  			rstream.packed(rngstream) <- TRUE
      
  		component_fits <- parallel_lapply_method(seq_len(crf_control$M_bag),
  			function(ind, data_concat, N, sequence_bag_groups, crf_control, rngstream, rstream_offset, resume_from) {
  				crf_control$save_file <- paste0(crf_control$save_path, crf_control$save_filename_base, "_baggroup", ind, ".Rdata")
  				validate_subjects <- which(!(seq_len(N) %in% sequence_bag_groups[[ind]]))
  				lccrf_one_validation_set(data_concat = data_concat, beta_penalty_factor = 0,
  					train_subjects = sequence_bag_groups[[ind]], validate_subjects = validate_subjects, crf_control = crf_control,
  					rngstream = rngstream, rstream_offset = rstream_offset + ind, resume_from = resume_from)
  			},
  			data_concat = data_concat, N = N, sequence_bag_groups = sequence_bag_groups, crf_control = crf_control, rngstream = rngstream, rstream_offset = 0, resume_from = resume_from)
  
  		if(identical(crf_control$rng_method, "rstream")) {
  			rstream.packed(rngstream) <- FALSE
			# would need to advance rngstream by crf_control$M_bag indices here if more rng were being performed below -- but it is not.
		}

  		temp_result <- list(component_fits = component_fits, sequence_bag_groups = sequence_bag_groups,
			oob_prop_correct = NA, oob_prop_correct_by_seq = NA)
  	} else if(identical(crf_control$bag_method, "none")) {
  		if(identical(crf_control$crossval_method, "fit-one")) {
  			crf_control$save_file <- paste0(crf_control$save_path, crf_control$save_filename_base, "_valgroup1.Rdata")
      
  			if(identical(crf_control$rng_method, "rstream")) {
  			  rstream.packed(rngstream) <- TRUE
  			}
      
  			train_subjects <- which(!(seq_len(N) %in% crossval_groups[[1]]))
  			component_fits <- list(lccrf_one_validation_set(data_concat = data_concat, beta_penalty_factor = 0,
  				train_subjects = train_subjects, validate_subjects = crossval_groups[[1]], crf_control = crf_control, rngstream = rngstream, rstream_offset = 0, resume_from = resume_from))
  
  			if(identical(crf_control$rng_method, "rstream"))
  			  rstream.packed(rngstream) <- FALSE
			# would need to advance rngstream by 1 index here if more rng were being performed below -- but it is not.
  
  			temp_result <- list(component_fits = component_fits, crossval_groups = crossval_groups[1])
  		} else if(identical(crf_control$crossval_method, "fit-all")) {
  			if(identical(crf_control$rng_method, "rstream"))
  			  rstream.packed(rngstream) <- TRUE
      
  			component_fits <- parallel_lapply_method(crossval_group_inds,
  				function(ind, data_concat, N, crossval_groups, crf_control, rngstream, rstream_offset, resume_from) {
  					crf_control$save_file <- paste0(crf_control$save_path, crf_control$save_filename_base, "_valgroup", ind, ".Rdata")
  					train_subjects <- which(!(seq_len(N) %in% crossval_groups[[ind]]))
  					lccrf_one_validation_set(data_concat = data_concat, beta_penalty_factor = 0,
  						train_subjects = train_subjects, validate_subjects = crossval_groups[[ind]],
  						crf_control = crf_control, rngstream = rngstream, rstream_offset = rstream_offset + ind, resume_from = resume_from)
  				},
  				data_concat = data_concat, N = N, crossval_groups = crossval_groups, crf_control = crf_control, rngstream = rngstream, rstream_offset = 0, resume_from = resume_from)
  
  			if(identical(crf_control$rng_method, "rstream"))
  			  rstream.packed(rngstream) <- FALSE
			# would need to advance rngstream by crf_control$M_bag indices here if more rng were being performed below -- but it is not.
  
  			temp_result <- list(component_fits = component_fits, crossval_groups = crossval_groups)
  		} else {
  			stop("Invalid crossval_method: must be fit-one or fit-all")
  		}
  	} else {
  		stop("Invalid bag_method: must be none or sequence")
  	} # end if bag_method
  } else if(identical(crf_control$fit_method, "parametric-L2-penalized-MLE")) {
    if(identical(crf_control$beta_penalty_factor, "crossval-select")) {
  		if(identical(crf_control$rng_method, "rstream"))
  			rstream.packed(rngstream) <- TRUE

    	crf_control$beta_penalty_factor <- list(crossval_select_beta_penalty_factor(data_concat, prev_log_obs_probs_by_state, first_term, prev_iter_log_lik, rngstream, rstream_offset = 0, crf_control, crossval_groups),
                                           "crossval-select")

  		if(identical(crf_control$rng_method, "rstream")) {
  			rstream.packed(rngstream) <- FALSE
			for(i in seq_len(crf_control$K_crossval)) {
				rstream.nextsubstream(rngstream)
			}
		}
  	} else if(!is.numeric(crf_control$beta_penalty_factor)) {
      stop("Invalid value for crf_control$beta_penalty_parameter: must be either numeric or \"crossval-select\"")
  	}
  	
  	temp_result <- list(component_fits = list(lccrf_one_validation_set(data_concat = data_concat, beta_penalty_factor = crf_control$beta_penalty_factor[[1]],
		  	train_subjects = seq_along(unique(data_concat$subject)), validate_subjects = c(), crf_control = crf_control, rngstream = rngstream, rstream_offset = 0,
			  resume_from = resume_from)))
  }

	if(crf_control$fit_method %in%  c("parametric-boost", "parametric-L2-penalized-MLE")) {
		temp_result$cent <- cent
		temp_result$scal <- scal
	}
  
	temp_result$crf_control <- crf_control
	temp_result$orig_y_levels <- orig_y_levels

	return(temp_result)
}


lccrf_one_validation_set <- function(data_concat, beta_penalty_factor, train_subjects, validate_subjects, crf_control, rngstream, rstream_offset, resume_from) {
	if(!is.null(resume_from))
		stop("resume functionality is not yet implemented")
  
	# handle rng
	if(identical(crf_control$rng_method, "rstream")) {
		rstream.packed(rngstream) <- FALSE
    
		# advance to substream specific to this crossval group
		for(i in seq_len(rstream_offset))
			rstream.nextsubstream(rngstream)
	  
		# set rstream as the RNG used by R, using the rngstream object
		rstream.RNG(rngstream)
	}
	
	# number of subjects
	N_validate <- length(validate_subjects)
	N_train <- length(train_subjects)
	N <- N_validate + N_train

	# number of time points for each subject and in total
	T_sub <- as.vector(table(data_concat$subject))
	if(N_validate > 0) {
		T_sub_validate <- T_sub[validate_subjects]
	} else {
    T_sub_validate <- 0
	}
	T_sub_train <- T_sub[train_subjects]

	T_total <- sum(T_sub_validate) + sum(T_sub_train)

	# number of variables
	D <- ncol(data_concat$X)

	# number of states
	S <- length(levels(data_concat$y))

	# if specified, add X^2 terms
	if(crf_control$quadratic)
		data_concat$X <- cbind(data_concat$X, (data_concat$X[, 1 + seq_len(D - 1)])^2)

    # estimate log pi based on full data set
    log_marginal_class_probs <- sapply(levels(data_concat$y), function(lvl) {log(sum(data_concat$y == lvl)) - log(T_total)})
    
	# get a C representation of a CRF
	CRF_Xptr <- initialize_new_CRF(D, S, log_marginal_class_probs)


	# separate data into training and validation subsets
	train_inds <- unlist(lapply(train_subjects, function(subj) which(data_concat$subject == subj)))
	
	train_data_concat <- list(X = data_concat$X[train_inds, , drop = FALSE],
		y = data_concat$y[train_inds],
		subject = as.integer(rep(seq_len(N_train), T_sub_train)))

	train_data_split <- lapply(unique(train_data_concat$subject), function(subj) {
		list(y = train_data_concat$y[train_data_concat$subject == subj],
			X = train_data_concat$X[train_data_concat$subject == subj, , drop = FALSE])
	})

	if(N_validate > 0) {
		validate_inds <- data_concat$subject %in% validate_subjects

		validate_data_concat <- list(X = data_concat$X[validate_inds, , drop = FALSE],
			y = data_concat$y[validate_inds],
			subject = as.integer(rep(seq_len(N_validate), T_sub_validate)))

		validate_data_split <- lapply(unique(validate_data_concat$subject), function(subj) {
			list(y = validate_data_concat$y[validate_data_concat$subject == subj],
				X = validate_data_concat$X[validate_data_concat$subject == subj, , drop = FALSE])
		})
	}

	# variables to track observation probabilities
	train_prev_log_obs_probs_by_state <- lapply(train_data_split, function(obs_dat) {
		matrix(0, nrow = nrow(obs_dat$X), ncol = S)
	})
  
	if(N_validate > 0) {
	  validate_prev_log_obs_probs_by_state <- lapply(validate_data_split, function(obs_dat) {
		  matrix(0, nrow = nrow(obs_dat$X), ncol = S)
  	})
	}

	if(is.null(resume_from)) {
		m <- 1L

    if(N_validate > 0 && !is.null(crf_control$M_boost)) {
  		validate_log_lik_trace <- rep(NA, crf_control$M_boost)
  		validate_prop_correct_trace <- rep(NA, crf_control$M_boost)
  		train_log_lik_trace <- rep(NA, crf_control$M_boost)
    } else {
      validate_log_lik_trace <- NA
      validate_prop_correct_trace <- NA
      train_log_lik_trace <- NA
    }
    
		prev_iter_log_lik <- -Inf
	}

	# estimate the parameters for the additive components
	prev_omega_hat <- NULL
	omega_hat <- list()
	active_vars <- list()
	log_alpha_hat <- list()
	active_cols <- list()
	beta_hat <- list()

	save_count <- 0L

	while(m <= crf_control$M_boost
			|| ifelse(length(which.max(validate_prop_correct_trace)) > 0,
					length(validate_prop_correct_trace) - which.max(validate_prop_correct_trace) <= crf_control$M_boost_search_threshold,
					FALSE)) {
		print(paste("m = ", m, sep = ""))

		next_component <- lccrf_est_one_component_parametric(CRF_Xptr = CRF_Xptr, data_concat = train_data_concat, data_split = train_data_split,
			D = D, prev_omega_hat = prev_omega_hat, beta_penalty_factor = beta_penalty_factor, prev_log_obs_probs_by_state = train_prev_log_obs_probs_by_state,
			first_term = identical(m, 1L), prev_iter_log_lik = prev_iter_log_lik, crf_control = crf_control)
		
		if(!is.null(next_component)) {
			prev_omega_hat <- next_component$omega_hat
			omega_hat[[m]] <- next_component$omega_hat
			active_vars[[m]] <- next_component$active_vars
			log_alpha_hat[[m]] <- next_component$log_alpha_hat
      
			active_cols[[m]] <- next_component$active_cols
			beta_hat[[m]] <- next_component$beta_hat

			prev_iter_log_lik <- next_component$log_lik
			train_log_lik_trace[m] <- next_component$log_lik
			train_prev_log_obs_probs_by_state <- next_component$combined_log_obs_probs_by_state

			if(N_validate > 0) {
				# omegas in proper format for lccrf_update_log_lik
				omegas <- calc_omegas_from_omega_hat(omega_hat[[m]], S, crf_control)

				validate_result_update <- lccrf_update_log_lik_parametric(betas = beta_hat[[m]], omegas = omegas,
 					active_cols = which(active_cols[[m]]), beta_penalty_factor = beta_penalty_factor, CRF_Xptr = CRF_Xptr,
 					data_split = validate_data_split, log_alpha = log_alpha_hat[[m]],
 					prev_log_obs_probs_by_state = validate_prev_log_obs_probs_by_state, crf_control = crf_control, retall = 1L)

				validate_log_lik_trace[m] <- validate_result_update[[1]]

				validate_prev_log_obs_probs_by_state <- validate_result_update[[2]]
      
				y_pred <- predict_lccrf_given_obs_probs_marginal(CRF_Xptr = CRF_Xptr, observedHMMs = validate_data_split,
					log_obs_probs_by_state = validate_prev_log_obs_probs_by_state, return_class_probs = FALSE)
				y_pred <- rbind.fill.matrix(y_pred)
      
				num_correct <- sum(y_pred == validate_data_concat$y)
				validate_prop_correct_trace[m] <- num_correct / length(validate_data_concat$y)

#				cat(paste0("validation log lik is ", validate_log_lik_trace[m]))
			}
		} else {
			m <- crf_control$M_boost + 1
			print(m)
			if(length(which.max(validate_prop_correct_trace)) > 0)
				validate_prop_correct_trace <- c(validate_prop_correct_trace, rep(NA, crf_control$M_boost_search_threshold + 1))
		}
    
		m <- m + 1L
	}

    if(N_validate > 0) {
  		return(list(log_pi_hat = log_marginal_class_probs, omega_hat = omega_hat, beta_hat = beta_hat, active_vars = active_vars, active_cols = active_cols, log_alpha_hat = log_alpha_hat, validate_subjects = validate_subjects, validate_log_lik_trace = validate_log_lik_trace, validate_prop_correct_trace = validate_prop_correct_trace, train_log_lik_trace = train_log_lik_trace))
	} else {
	    return(list(log_pi_hat = log_marginal_class_probs, omega_hat = omega_hat, beta_hat = beta_hat, active_vars = active_vars, active_cols = active_cols, log_alpha_hat = log_alpha_hat, validate_subjects = validate_subjects, train_log_lik_trace = train_log_lik_trace))
	}
}


lccrf_est_one_component_parametric <- function(CRF_Xptr, data_concat, data_split, D, prev_omega_hat, beta_penalty_factor, prev_log_obs_probs_by_state, first_term, prev_iter_log_lik, crf_control) {
	# it is assumed that the data have been augmented so that the first column of X is a vector of 1's.

	# number of time points for each subject and in total
	T_sub <- as.vector(table(data_concat$subject))
	T_total <- length(data_concat$subject)

	# number of states
	S <- ncol(prev_log_obs_probs_by_state[[1]])
  
	full_D <- ncol(data_concat$X)

	# no shrinkage
	log_alpha_hat <- 0

	# estimate new component parameters	
	log_lik <- prev_iter_log_lik
	attempt_num <- 0
	vars_not_tried <- seq_len(D - 1)

	log_marginal_class_probs <- sapply(levels(data_concat$y), function(lvl) {log(sum(data_concat$y == lvl)) - log(T_total)})

	max_vars <- as.integer(min(crf_control$num_active_vars[[1]], D - 1))

	while(attempt_num < crf_control$max_attempts && ifelse(identical(max_vars, 1L), length(vars_not_tried) > 1, TRUE) && prev_iter_log_lik >= log_lik) {
		cat('.')
		attempt_num <- attempt_num + 1
		
		if(first_term) {
		  #			data_concat <- list(y = rbind.fill.matrix(lapply(data_split, function(comp) as.integer(comp$y))))
		  #			data_concat$subject <- rep(seq_along(data_split), sapply(data_split, function(comp) length(comp$y)))
		  if(crf_control$reduced_trans_mat_parameterization) {
			  log_trans_mat <- estimate_transition_matrix_reduced_parameterization_states_observed(data_concat$y, data_concat$subject, S, log = TRUE)
			  omega_hat <- log_trans_mat[1, 1]
		  } else {
			  log_trans_mat <- estimate_transition_matrix_full_parameterization_states_observed(data_concat$y, data_concat$subject, S, log = TRUE)
			  omega_hat <- as.vector(t(log_trans_mat))
			  omega_hat[omega_hat < -100] <- -100
			  omega_hat <- omega_hat - omega_hat[length(omega_hat)]
			  omega_hat <- omega_hat[seq_len(length(omega_hat) - 1)]
			  omega_hat[omega_hat < -100] <- -100
		  }
		} else {
		  omega_hat <- prev_omega_hat
		}
		

		S_max <- S - 1

		active_vars <- rep(FALSE, D)

		if(identical(crf_control$active_var_search_method, "random")) {
			if(identical(max_vars, 1L)) {
				new_var_ind <- sample(length(vars_not_tried), 1)
				active_vars[c(1, 1 + vars_not_tried[new_var_ind])] <- TRUE
				vars_not_tried <- vars_not_tried[-new_var_ind]
			} else {
				active_vars[c(1, 1 + sample(D - 1, max_vars))] <- TRUE
			}
		} else if(identical(crf_control$active_var_search_method, "max-gradient")) {
#		  stop("The gradient-max active_var_search_method has not yet been implemented")
		  # calculate "combined gradient magnitude"

	    # omegas in proper format for lccrf_update_log_lik_new_betas_first
		  omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)
		  
		  beta_hat <- rep(0, S * full_D)
      
		  gr <- lccrf_update_gradient_neg_log_lik_new_betas_first_parametric(beta_hat[c(rep(TRUE, full_D * S_max), rep(FALSE, (S - S_max) * full_D))],
		                        omegas = omegas,
		                        active_cols = seq_len(full_D),
		                        beta_penalty_factor = beta_penalty_factor,
		                        CRF_Xptr = CRF_Xptr,
		                        data_split = data_split,
		                        S = S,
		                        S_max = S_max,
		                        D = D,
		                        log_alpha = log_alpha_hat,
		                        prev_log_obs_probs_by_state = prev_log_obs_probs_by_state,
		                        crf_control = crf_control,
		                        diff_args = list(first_term = first_term, beta_diff_cols = seq_len(full_D), diff_wrt_omegas = FALSE))
		                        
		  cgm <- apply(abs(matrix(gr, nrow = S_max, byrow = TRUE)), 2, sum)[-1]
      
		  if(crf_control$quadratic) {
		    cgm <- apply(abs(matrix(cgm, nrow = 2 * S_max, byrow = TRUE)), 2, sum)
		  }
      
		  active_vars[c(1, 1 + which.max(cgm))] <- TRUE
		}

    if(crf_control$quadratic) {
      active_cols <- c(active_vars, active_vars[1 + seq_len(D - 1)])
    } else {
      active_cols <- active_vars
    }


		# initial values for betas
		beta_hat <- rep(0, S * full_D)
		if(identical(crf_control$beta_initialization, "random")) {
			# omegas in proper format for lccrf_update_log_lik_new_betas_first
			omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)
			
			# random start values
			beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))] <- rnorm(S_max * sum(active_cols))
		} else if(identical(crf_control$beta_initialization, "MLR")) {
			# get initial betas using multinomial logistic regression with offset
			# omegas in proper format
			omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)

			beta_hat[c(rep(active_cols, S - 1), rep(FALSE, 1 * full_D))] <- lccrf_compute_initial_betas_MLR(CRF_Xptr,
				active_cols = which(active_cols), data_concat = data_concat, data_split = data_split,
				T_sub = T_sub, S = S, omegas = omegas,
				prev_log_obs_probs_by_state = prev_log_obs_probs_by_state, crf_control = crf_control)
		}

		# optimize betas only
		optim_result <- optim(par = beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))],
													fn = lccrf_update_neg_log_lik_new_betas_first_parametric,
													gr = lccrf_update_gradient_neg_log_lik_new_betas_first_parametric,
													#				gr = NULL,
													omegas = omegas,
													active_cols = which(active_cols),
													beta_penalty_factor = beta_penalty_factor,
													CRF_Xptr = CRF_Xptr,
													data_split = data_split,
													S = S,
													S_max = S_max,
													D = D,
													log_alpha = log_alpha_hat,
													prev_log_obs_probs_by_state = prev_log_obs_probs_by_state,
													crf_control = crf_control,
													diff_args = list(first_term = first_term, beta_diff_cols = which(active_cols), diff_wrt_omegas = FALSE),
													method = crf_control$optim_method,
													#		lower = -10000,
													#		upper = 10000,
													control = list(),
													hessian = FALSE)
		
		beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))] <- optim_result$par


		
      if(crf_control$update_transition_matrix) { # optimize all parameters including transition matrix
  		optim_result <- optim(par = c(omega_hat, beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))]),
  			fn = lccrf_update_neg_log_lik_new_params_first_parametric,
  			gr = lccrf_update_gradient_neg_log_lik_new_params_first_parametric,
  #			gr = NULL,
  			active_cols = which(active_cols),
  			beta_penalty_factor = beta_penalty_factor,
  			CRF_Xptr = CRF_Xptr,
  			data_split = data_split,
  			S = S,
  			S_max = S_max,
  			D = D,
  			log_alpha = log_alpha_hat,
  			prev_log_obs_probs_by_state = prev_log_obs_probs_by_state,
  			crf_control = crf_control,
  			diff_args = list(first_term = first_term, beta_diff_cols = which(active_cols), diff_wrt_omegas = TRUE),
  			method = crf_control$optim_method,
  	#		lower = -10000,
  	#		upper = 10000,
  			control = list(),
  			hessian = FALSE)
  		
  		if(crf_control$reduced_trans_mat_parameterization) {
  			omega_hat <- optim_result$par[1]
  			beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))] <- optim_result$par[seq(from = 2, length = S_max * sum(active_cols), by = 1)]
  		} else {
  			omega_hat <- optim_result$par[seq_len(S^2 - 1)]
  			beta_hat[c(rep(active_cols, S_max), rep(FALSE, (S - S_max) * full_D))] <- optim_result$par[seq(from = S^2, length = S_max * sum(active_cols), by = 1)]
  		}
      }

		log_lik <- -1 * optim_result$value
	}

	# return null if we did not find an improvement
	if(log_lik <= prev_iter_log_lik)
		return(NULL)

	# update parameters, log_lik and combined_log_obs_probs_by_state -- this step may not be necessary now that shrinkage has been removed
	
	# omegas in proper format for lccrf_update_log_lik
	omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)

	temp <- lccrf_update_log_lik_parametric(betas = beta_hat, omegas = omegas,
		active_cols = which(active_cols), beta_penalty_factor = beta_penalty_factor, CRF_Xptr = CRF_Xptr, data_split = data_split,
		log_alpha = log_alpha_hat, prev_log_obs_probs_by_state = prev_log_obs_probs_by_state, crf_control = crf_control, retall = 1L)

	log_lik <- temp[[1]]
	combined_log_obs_probs_by_state <- temp[[2]]

	return(list(omega_hat = omega_hat, active_vars = active_vars, active_cols = active_cols, beta_hat = beta_hat, log_alpha_hat = log_alpha_hat, log_lik = log_lik, combined_log_obs_probs_by_state = combined_log_obs_probs_by_state))
}

crossval_select_beta_penalty_factor <- function(data_concat, prev_log_obs_probs_by_state, first_term, prev_iter_log_lik, rngstream, rstream_offset, crf_control, crossval_groups) {
	data_split <- lapply(unique(data_concat$subject), function(subj) {
		list(y = data_concat$y[data_concat$subject == subj],
			X = data_concat$X[data_concat$subject == subj, , drop = FALSE])
	})

	# optimize beta_penalty_factor
  optim_result <- optim(par = 1,
#                        fn = lccrf_crossval_est_neg_prop_correct_by_beta_penalty_factor,
                        fn = lccrf_crossval_est_neg_val_log_lik_by_beta_penalty_factor,
                        gr = NULL,
                        data_concat = data_concat,
          						  data_split = data_split,
          						  crossval_group_inds = seq_along(crossval_groups),
          						  crossval_groups = crossval_groups,
          						  rngstream = rngstream,
          						  rstream_offset = rstream_offset,
          						  crf_control = crf_control,
#                        method = "Brent",
												method = "L-BFGS-B",
												lower = 0,
#                        upper = 10000,
                        control = list(),
                        hessian = FALSE)
#  browser()
	return(optim_result$par)
}

lccrf_crossval_est_neg_prop_correct_by_beta_penalty_factor <- function(beta_penalty_factor, data_concat, data_split, crossval_group_inds, crossval_groups, rngstream, rstream_offset, crf_control) {
  if(identical(crf_control$parallel_method, "none")) {
    crossval_lapply_method <- lapply
  } else if(identical(crf_control$parallel_method, "snowfall")) {
    crossval_lapply_method <- sfLapply
  } else {
    stop("Invalid parallel_method: must be none or snowfall")
  }

  D <- ncol(data_concat$X)
  
  rstream.packed(rngstream) <- TRUE

  crossval_nums_correct <- crossval_lapply_method(crossval_group_inds,
    function(ind, data_concat, crossval_groups, beta_penalty_factor, rngstream, rstream_offset, crf_control) {
  	  # get reduced data not including crossval_groups[ind], then pass in train_subjects and validate_subjects based on M_selection_crossval_subgroups
  	  crossval_group_subjs <- unlist(crossval_groups[-ind])
  	  crossval_group_obs_inds <- unlist(lapply(crossval_group_subjs, function(subj) which(data_concat$subject == subj)))
  	  crossval_group_T_sub <- sapply(crossval_group_subjs, function(subj) sum(data_concat$subject == subj))
	
  	  crossval_group_data_concat <- list(X = data_concat$X[crossval_group_obs_inds, , drop = FALSE],
  		y = data_concat$y[crossval_group_obs_inds],
  		subject = as.integer(rep(seq_along(crossval_group_subjs), crossval_group_T_sub)))

  	  crossval_train_subjs <- seq_along(crossval_group_subjs)
      crf_control$save_freq <- Inf
      component_fit <- lccrf_one_validation_set(data_concat = crossval_group_data_concat, beta_penalty_factor = beta_penalty_factor,
        train_subjects = crossval_train_subjs, validate_subjects = c(),
        rngstream = rngstream, rstream_offset = rstream_offset + ind,
        crf_control = crf_control, resume_from = NULL)

      num_correct <- sum(sapply(seq_along(crossval_groups[[ind]]), function(group_sub) {
        y_pred_one_sub <- predict_lccrf(data_split[crossval_groups[[ind]][group_sub]],
                                    fit = list(crf_control = crf_control, component_fits = list(component_fit), cent = rep(0, D), scal = rep(1, D), orig_y_levels = levels(data_split[[1]]$y)))[[1]]
        sum(data_split[[crossval_groups[[ind]][group_sub]]]$y == y_pred_one_sub)
      }))

      return(num_correct)
    },
    data_concat = data_concat, crossval_groups = crossval_groups, beta_penalty_factor = beta_penalty_factor, rngstream = rngstream, rstream_offset = rstream_offset, crf_control = crf_control)

  crossval_nums_correct <- unlist(crossval_nums_correct)

  rstream.packed(rngstream) <- FALSE

  # IF ANYTHING BELOW HERE REQUIRED RNG, I MIGHT HAVE TO RESET rngstream -- BUT IT DOESN'T
  
  return(-1 * sum(crossval_nums_correct) / length(data_concat$y))
}


lccrf_crossval_est_neg_val_log_lik_by_beta_penalty_factor <- function(beta_penalty_factor, data_concat, data_split, crossval_group_inds, crossval_groups, rngstream, rstream_offset, crf_control) {
  if(identical(crf_control$parallel_method, "none")) {
    crossval_lapply_method <- lapply
  } else if(identical(crf_control$parallel_method, "snowfall")) {
    crossval_lapply_method <- sfLapply
  } else {
    stop("Invalid parallel_method: must be none or snowfall")
  }

  D <- ncol(data_concat$X)
  
  rstream.packed(rngstream) <- TRUE

  crossval_log_liks <- crossval_lapply_method(crossval_group_inds,
    function(ind, data_concat, crossval_groups, beta_penalty_factor, rngstream, rstream_offset, crf_control) {
  	  # get reduced data not including crossval_groups[ind], then pass in train_subjects and validate_subjects based on M_selection_crossval_subgroups
  	  crossval_group_subjs <- unlist(crossval_groups[-ind])
  	  crossval_group_obs_inds <- unlist(lapply(crossval_group_subjs, function(subj) which(data_concat$subject == subj)))
  	  crossval_group_T_sub <- sapply(crossval_group_subjs, function(subj) sum(data_concat$subject == subj))
	
  	  crossval_group_data_concat <- list(X = data_concat$X[crossval_group_obs_inds, , drop = FALSE],
  		y = data_concat$y[crossval_group_obs_inds],
  		subject = as.integer(rep(seq_along(crossval_group_subjs), crossval_group_T_sub)))

  	  crossval_train_subjs <- seq_along(crossval_group_subjs)
      crf_control$save_freq <- Inf
      component_fit <- lccrf_one_validation_set(data_concat = crossval_group_data_concat, beta_penalty_factor = beta_penalty_factor,
        train_subjects = crossval_train_subjs, validate_subjects = c(),
        rngstream = rngstream, rstream_offset = rstream_offset + ind,
        crf_control = crf_control, resume_from = NULL)
    
      log_lik <- sum(sapply(seq_along(crossval_groups[[ind]]), function(group_sub) {
        S <- length(levels(data_split[[1]]$y))
        T_sub <- sapply(data_split, function(comp) length(comp$y))
        T_total <- sum(T_sub)
        
        D <- ncol(data_split[[1]]$X)
        
        # get a C representation of a CRF
        CRF_Xptr <- initialize_new_CRF(D, S)
        
        # combined_log_obs_probs
        prev_log_obs_probs_by_state <- lapply(data_split[crossval_groups[[ind]][group_sub]], function(obs_dat) {
          matrix(0, nrow = nrow(obs_dat$X), ncol = S)
        })
    
        omegas <- calc_omegas_from_omega_hat(component_fit$omega_hat[[1]], S, crf_control)
        
        return(lccrf_update_log_lik_parametric(betas = component_fit$beta_hat[[1]], omegas = omegas,
          active_cols = component_fit$active_cols[[1]], beta_penalty_factor = 0, CRF_Xptr = CRF_Xptr,
          data_split = data_split[crossval_groups[[ind]][group_sub]], log_alpha = 0,
          prev_log_obs_probs_by_state = prev_log_obs_probs_by_state, crf_control = crf_control, retall = FALSE))
      }))

      return(log_lik)
    },
    data_concat = data_concat, crossval_groups = crossval_groups, beta_penalty_factor = beta_penalty_factor, rngstream = rngstream, rstream_offset = rstream_offset, crf_control = crf_control)

  crossval_log_liks <- unlist(crossval_log_liks)

  rstream.packed(rngstream) <- FALSE

  # IF ANYTHING BELOW HERE REQUIRED RNG, I MIGHT HAVE TO RESET rngstream -- BUT IT DOESN'T
  cat(-1 * sum(crossval_log_liks))
  return(-1 * sum(crossval_log_liks))
}


## likelihood functions

lccrf_parse_params_vector <- function(params, active_cols, S, S_max, D, full_D, crf_control) {
	if(crf_control$reduced_trans_mat_parameterization) {
		num_omegas <- 1
	} else {
		num_omegas <- S^2 - 1
	}
	num_betas_per_state <- length(active_cols)

	# get omegas
	omega_hat <- params[seq_len(num_omegas)]
	omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)

	# get betas
	betas <- rep(0, S * full_D)
	for(s in seq_len(S_max))
		betas[(s - 1) * full_D + active_cols] <- params[num_omegas + (s - 1) * num_betas_per_state + seq_len(num_betas_per_state)]

	return(list(omegas = omegas, betas = betas))
}

lccrf_update_log_lik_new_params_first_parametric <- function(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S, S_max,
		D, log_alpha, prev_log_obs_probs_by_state, crf_control, retall, diff_args) {

	full_D <- ncol(data_split[[1]]$X)

	parsed_params <- lccrf_parse_params_vector(params, active_cols, S, S_max, D, full_D, crf_control)

	lccrf_update_log_lik_parametric(parsed_params$betas, parsed_params$omegas,
		active_cols, beta_penalty_factor = beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
		prev_log_obs_probs_by_state, crf_control, retall)
}

lccrf_update_neg_log_lik_new_params_first_parametric <- function(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S,
		S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, diff_args) {
	return(-1 * lccrf_update_log_lik_new_params_first_parametric(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S, S_max, D, log_alpha,
		prev_log_obs_probs_by_state, crf_control, retall = 0L, diff_args))
}


lccrf_update_log_lik_new_betas_first_parametric <- function(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
	S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, retall, diff_args) {
	full_D <- ncol(data_split[[1]]$X)

	num_betas_per_state <- length(active_cols)

	if(length(betas_vec) != S_max * num_betas_per_state)
		stop("params must have length S^2 - 1 + (S - 1) * length(active_cols) + D")

	# get betas
	betas <- rep(0, S * full_D)
	for(s in seq_len(S_max))
		betas[(s - 1) * full_D + active_cols] <- betas_vec[(s - 1) * num_betas_per_state + seq_len(num_betas_per_state)]

	lccrf_update_log_lik_parametric(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
		prev_log_obs_probs_by_state, crf_control, retall)
}

lccrf_update_neg_log_lik_new_betas_first_parametric <- function(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
	S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, diff_args) {
	return(-1 * lccrf_update_log_lik_new_betas_first_parametric(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr,
		data_split, S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, retall = 0L, diff_args))
}


lccrf_update_log_lik_parametric <- function(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
	prev_log_obs_probs_by_state, crf_control, retall) {

	result <- globalSeqCRF_update_log_lik(betas, omegas, active_cols, CRF_Xptr, data_split,
		prev_log_obs_probs_by_state, retall)

	if(identical(crf_control$fit_method, "parametric-L2-penalized-MLE")) {
    if(retall) {
      result[[1]] <- result[[1]] - beta_penalty_factor * sum(betas^2)
    } else {
  		result <- result - beta_penalty_factor * sum(betas^2)
    }
	}

	return(result)
}

lccrf_update_neg_log_lik_parametric <- function(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
	prev_log_obs_probs_by_state, crf_control) {
	return(-1 * lccrf_update_log_lik_parametric(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
		log_alpha, prev_log_obs_probs_by_state, crf_control, retall = 0L))
}

lccrf_update_log_lik_given_log_obs_probs_omega_hat_first_parametric <- function(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
    crf_control, diff_args, time_scale, retall) {
  # get omegas
  omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)
  set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
  
  CRF_log_lik_given_obs_probs(CRF_Xptr, data_split, log_obs_probs_by_state, time_scale, retall)
}

lccrf_update_neg_log_lik_given_log_obs_probs_omega_hat_first_parametric <- function(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
                                                                                       crf_control, diff_args, time_scale) {
  return(-1 * lccrf_update_log_lik_given_log_obs_probs_omega_hat_first_parametric(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
                                                                                         crf_control, diff_args, time_scale, retall = 0L))
}

CRF_log_lik_given_obs_probs <- function(CRF_Xptr, observedHMMs, log_obs_probs_by_state, time_scale, retall) {
	if(identical(time_scale, "marginal")) {
		return(.Call("marginal_CRF_log_lik_given_obs_probs_C", CRF_Xptr, observedHMMs, log_obs_probs_by_state, as.integer(retall)))
	} else if(identical(time_scale, "sequence")) {
		return(.Call("CRF_log_lik_given_obs_probs_C", CRF_Xptr, observedHMMs, log_obs_probs_by_state))
	} else {
		stop("time_scale must be either \"marginal\" or \"sequence\"")
	}
}

CRF_neg_log_lik_given_obs_probs <- function(CRF_Xptr, observedHMMs, log_obs_probs_by_state, time_scale, retall) {
	return(-1 * CRF_log_lik_given_obs_probs(CRF_Xptr, observedHMMs, log_obs_probs_by_state, time_scale, retall))
}



## gradient functions
lccrf_update_gradient_log_lik_new_params_first_parametric <- function(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S, S_max, D, log_alpha,
	prev_log_obs_probs_by_state, crf_control, retall, diff_args) {

	full_D <- ncol(data_split[[1]]$X)

	parsed_params <- lccrf_parse_params_vector(params, active_cols, S, S_max, D, full_D, crf_control)

	lccrf_update_gradient_log_lik_parametric(parsed_params$betas, parsed_params$omegas,
		active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
		prev_log_obs_probs_by_state, crf_control, diff_args)
}

lccrf_update_gradient_neg_log_lik_new_params_first_parametric <- function(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S, S_max, D, log_alpha,
	prev_log_obs_probs_by_state, crf_control, diff_args) {
	return(-1 * lccrf_update_gradient_log_lik_new_params_first_parametric(params, active_cols, beta_penalty_factor, CRF_Xptr, data_split, S, S_max, D, log_alpha,
		prev_log_obs_probs_by_state, crf_control, retall = 0L, diff_args))
}


lccrf_update_gradient_log_lik_new_betas_first_parametric <- function(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
	S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, retall, diff_args) {
	full_D <- ncol(data_split[[1]]$X)

	num_betas_per_state <- length(active_cols)

	if(length(betas_vec) != S_max * num_betas_per_state)
		stop("params must have length S^2 - 1 + (S - 1) * length(active_cols) + D")

	# get betas
	betas <- rep(0, S * full_D)
	for(s in seq_len(S_max))
		betas[(s - 1) * full_D + active_cols] <- betas_vec[(s - 1) * num_betas_per_state + seq_len(num_betas_per_state)]

	lccrf_update_gradient_log_lik_parametric(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
		prev_log_obs_probs_by_state, crf_control, diff_args)
}

lccrf_update_gradient_neg_log_lik_new_betas_first_parametric <- function(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
	S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, diff_args) {
	return(-1 * lccrf_update_gradient_log_lik_new_betas_first_parametric(betas_vec, omegas, active_cols, beta_penalty_factor, CRF_Xptr,
		data_split, S, S_max, D, log_alpha, prev_log_obs_probs_by_state, crf_control, retall = 0L, diff_args))
}


lccrf_update_gradient_log_lik_parametric <- function(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
	prev_log_obs_probs_by_state, crf_control, diff_args) {

	result <- globalSeqCRF_update_gradient_log_lik(betas, omegas, active_cols, CRF_Xptr, data_split,
		prev_log_obs_probs_by_state, diff_args$beta_diff_cols, diff_args$diff_wrt_omegas, as.integer(crf_control$reduced_trans_mat_parameterization))

	if(identical(crf_control$fit_method, "parametric-L2-penalized-MLE")) {
		S <- length(omegas)
		num_betas_per_state <- length(betas) / S
		penalty_terms <- beta_penalty_factor * 2 * betas[rep(diff_args$beta_diff_cols, S - 1) + 
			rep(seq_len(S - 1) - 1, each = length(diff_args$beta_diff_cols)) * num_betas_per_state]
		
		if(diff_args$diff_wrt_omegas) {
			if(crf_control$reduced_trans_mat_parameterization) {
				num_omegas <- 1
			} else {
				num_omegas <- S^2 - 1
			}
		} else {
			num_omegas <- 0
		}
		
		beta_inds <- num_omegas + seq_len((S - 1) * length(diff_args$beta_diff_cols))
		result[beta_inds] <- result[beta_inds] - penalty_terms
	}

	return(result)
}

lccrf_update_gradient_neg_log_lik_parametric <- function(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split, log_alpha,
	prev_log_obs_probs_by_state, crf_control) {
	return(-1 * lccrf_update_log_lik_parametric(betas, omegas, active_cols, beta_penalty_factor, CRF_Xptr, data_split,
		log_alpha, prev_log_obs_probs_by_state, crf_control, retall = 0L))
}

lccrf_update_gradient_log_lik_given_log_obs_probs_omega_hat_first_parametric <- function(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
                                                                                   crf_control, diff_args, time_scale) {
  # get omegas
  omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)
  set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)

  if(identical(time_scale, "sequence")) {
    return(.Call("CRF_gradient_log_lik_wrt_omegas_given_obs_probs_C", CRF_Xptr, data_split, log_obs_probs_by_state, as.integer(crf_control$reduced_trans_mat_parameterization)))
  } else {
    stop("gradient wrt omega is not yet implemented for specified time scale")
  }
}

lccrf_update_gradient_neg_log_lik_given_log_obs_probs_omega_hat_first_parametric <- function(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
                                                                                       crf_control, diff_args, time_scale) {
  return(-1 * lccrf_update_gradient_log_lik_given_log_obs_probs_omega_hat_first_parametric(omega_hat, CRF_Xptr, data_split, S, log_obs_probs_by_state,
                                                                                     crf_control, diff_args, time_scale))
}


## prediction functions

predict_lccrf <- function(data_split, fit, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs = FALSE) {
	data_split <- lapply(data_split, function(comp) {
		comp$y <- factor(rep(1, nrow(comp$X)), levels = seq_along(fit$orig_y_levels))
		return(comp)
	})
	
	if(fit$crf_control$fit_method %in%  c("parametric-boost", "parametric-L2-penalized-MLE")) {
		predict_lccrf_parametric(data_split, fit, component_model_combine_method = component_model_combine_method, predict_method = predict_method, M_selection_method = M_selection_method, return_class_probs = return_class_probs)
	} else {
		stop("Invalid value for crf_control$fit_method.")
	}
}


predict_lccrf_parametric <- function(data_split, fit, component_model_combine_method = "equal-weight-lop", predict_method = "marginal", M_selection_method = "crossval-prop-correct", return_class_probs) {
  S <- length(levels(data_split[[1]]$y))
	T_sub <- sapply(data_split, function(comp) length(comp$y))
	T_total <- sum(T_sub)

	D <- ncol(data_split[[1]]$X)

	data_split <- lapply(data_split, function(comp) {
		if(!isTRUE(all.equal(comp$X[, 1], rep(1, nrow(comp$X)))) && fit$crf_control$add_intercept)
			comp$X <- cbind(rep(1, nrow(comp$X)), comp$X)

		comp$X <- scale(comp$X, center = fit$cent, scale = fit$scal)

		# if specified, add X^2 terms
		if((fit$crf_control)$quadratic)
			comp$X <- cbind(comp$X, (comp$X[, 1 + seq_len(ncol(comp$X) - 1)])^2)

		# add trailing column of 0's for L2 penalized method -- otherwise, coefficients would not be treated the same.
		if(identical(fit$crf_control$fit_method, "parametric-L2-penalized-MLE")) {
			comp$X <- cbind(comp$X, rep(0, nrow(comp$X)))
		}

		return(comp)
	})

	# get a C representation of a CRF
    CRF_Xptr <- initialize_new_CRF(D, S)

	# combined_log_class_probs
	combined_log_class_probs <- lapply(data_split, function(obs_dat) {
		matrix(0, nrow = nrow(obs_dat$X), ncol = S)
	})

	if(identical(component_model_combine_method, "equal-weight-marginalized-lop") && identical(predict_method, "marginal")) {
		num_component_fits <- length(fit$component_fits)
		
		for(fit_ind in seq_along(fit$component_fits)) {
			log_obs_probs_by_state <- lapply(data_split, function(obs_dat) {
				matrix(0, nrow = nrow(obs_dat$X), ncol = S)
			})

			# assume estimated log pi is the same in all components
			set_HMM_pi_from_log(CRF_Xptr, fit$component_fits[[fit_ind]]$log_pi_hat)

			if(length(fit$component_fits[[fit_ind]]$validate_log_lik_trace) <= 1) {
				M <- length(fit$component_fits[[fit_ind]]$beta_hat)
			} else {
				M <- which.max(fit$component_fits[[fit_ind]]$validate_log_lik_trace)
			}
      
			# get omegas
			omegas <- calc_omegas_from_omega_hat(omega_hat, S, crf_control)

			for(m in seq_len(M)) {
				temp <- lccrf_update_log_lik_parametric(betas = fit$component_fits[[fit_ind]]$beta_hat[[m]], omegas = omegas,
					active_cols = which(fit$component_fits[[fit_ind]]$active_cols[[m]]), beta_penalty_factor = 0, CRF_Xptr = CRF_Xptr,
					data_split = data_split, log_alpha = fit$component_fits[[fit_ind]]$log_alpha_hat[[m]],
					prev_log_obs_probs_by_state = log_obs_probs_by_state, crf_control = (fit$crf_control),
					retall = 1L)
				log_obs_probs_by_state <- temp[[2]]
			}

			new_log_class_probs <- CRF_log_lik_given_obs_probs(CRF_Xptr, data_split, log_obs_probs_by_state, time_scale = "marginal", retall = 1)[[2]]

			combined_log_class_probs <- lapply(seq_along(combined_log_class_probs), function(sub_ind) {
				combined_log_class_probs[[sub_ind]] + new_log_class_probs[[sub_ind]] / num_component_fits
			})
		}

#		combined_ps <- rbind.fill.matrix(combined_log_class_probs)
		if(return_class_probs) {
			return(list(predictions = lapply(combined_log_class_probs, function(combined_ps) {
				apply(combined_ps, 1, which.max)
			}),
			log_class_probs = combined_log_class_probs))
		} else {
			return(lapply(combined_log_class_probs, function(combined_ps) {
				apply(combined_ps, 1, which.max)
			}))
		}
	} else if(identical(component_model_combine_method, "equal-weight-lop") && identical(predict_method, "marginal")) {
		num_component_fits <- length(fit$component_fits)

		if(identical(M_selection_method, "crossval-lik")) {
    	cv_M <- sapply(fit$component_fits, function(cvfit) {
        M <- which.max(cvfit$validate_log_lik_trace)
        if(identical(integer(0), M))
          M <- length(cvfit$beta_hat)
        return(M)
      })
		} else if(identical(M_selection_method, "crossval-prop-correct")) {
		  cv_M <- sapply(fit$component_fits, function(cvfit) {
        M <- which.max(cvfit$validate_prop_correct_trace)
        if(identical(integer(0), M))
          M <- length(cvfit$beta_hat)
        return(M)
		  })
		} else if(identical(M_selection_method, "max")) {
      cv_M <- sapply(fit$component_fits, function(cvfit) {
        M <- length(cvfit$beta_hat)
        return(M)
      })
		} else {
      stop("Invalid M_selection_method: must be crossval-lik, crossval-prop-correct, or max")
		}

		# estimated log pi
		temp <- rbind.fill.matrix(lapply(fit$component_fits, function(cvfit) matrix(cvfit$log_pi_hat, nrow = 1)))
		combined_log_pi_hat <- apply(temp, 2, mean)
		set_HMM_pi_from_log(CRF_Xptr, combined_log_pi_hat)

		# estimated omegas
		temp <- rbind.fill.matrix(lapply(seq_len(num_component_fits), function(fit_ind) {
			matrix(fit$component_fits[[fit_ind]]$omega_hat[[cv_M[fit_ind]]], nrow = 1)
		}))
		combined_omega_hat <- apply(temp, 2, mean)
		omegas <- calc_omegas_from_omega_hat(combined_omega_hat, S, fit$crf_control)

		# estimated log "observation probabilities"
		log_obs_probs_by_state <- lapply(data_split, function(obs_dat) {
			matrix(0, nrow = nrow(obs_dat$X), ncol = S)
		})

		for(fit_ind in seq_along(fit$component_fits)) {
			for(m in seq_len(cv_M[fit_ind])) {
				temp <- lccrf_update_log_lik_parametric(betas = fit$component_fits[[fit_ind]]$beta_hat[[m]] / num_component_fits, omegas = omegas,
					active_cols = which(fit$component_fits[[fit_ind]]$active_cols[[m]]), beta_penalty_factor = 0, CRF_Xptr = CRF_Xptr,
					data_split = data_split, log_alpha = fit$component_fits[[fit_ind]]$log_alpha_hat[[m]],
					prev_log_obs_probs_by_state = log_obs_probs_by_state, crf_control = (fit$crf_control),
					retall = 1L)
				log_obs_probs_by_state <- temp[[2]]
			}
		}

		combined_log_class_probs <- CRF_log_lik_given_obs_probs(CRF_Xptr, data_split, log_obs_probs_by_state, time_scale = "marginal", retall = 1)[[2]]
#		combined_ps <- rbind.fill.matrix(combined_log_class_probs)
		if(return_class_probs) {
			return(list(predictions = lapply(combined_log_class_probs, function(combined_ps) {
				apply(combined_ps, 1, which.max)
			}),
			log_class_probs = combined_log_class_probs))
		} else {
			return(lapply(combined_log_class_probs, function(combined_ps) {
				apply(combined_ps, 1, which.max)
			}))
		}
	} else if(identical(component_model_combine_method, component_model_combine_method = "majority-vote-bag-lop-boost") && identical(predict_method, "marginal")) {
	  num_component_fits <- length(fit$component_fits)
	  
	  if(identical(M_selection_method, "crossval-lik")) {
	    cv_M <- sapply(fit$component_fits, function(cvfit) {
	      M <- which.max(cvfit$validate_log_lik_trace)
	      if(identical(integer(0), M))
	        M <- length(cvfit$beta_hat)
	      return(M)
	    })
	  } else if(identical(M_selection_method, "crossval-prop-correct")) {
	    cv_M <- sapply(fit$component_fits, function(cvfit) {
	      M <- which.max(cvfit$validate_prop_correct_trace)
	      if(identical(integer(0), M))
	        M <- length(cvfit$beta_hat)
	      return(M)
	    })
	  } else if(identical(M_selection_method, "max")) {
	    cv_M <- sapply(fit$component_fits, function(cvfit) {
	      M <- length(cvfit$beta_hat)
	      return(M)
	    })
	  } else {
	    stop("Invalid M_selection_method: must be crossval-lik, crossval-prop-correct, or max")
	  }
	  
	  # predicted state from each component model
	  y_pred_by_component_ind <- lapply(data_split, function(obs_dat) {
	    matrix(0, nrow = nrow(obs_dat$X), ncol = num_component_fits)
	  })
	  
	  for(fit_ind in seq_along(fit$component_fits)) {
	    # estimated log pi
	    set_HMM_pi_from_log(CRF_Xptr, fit$component_fits[[fit_ind]]$log_pi_hat)
	    
	    # estimated omegas
	    omegas <- calc_omegas_from_omega_hat(fit$component_fits[[fit_ind]]$omega_hat[[cv_M[fit_ind]]], S, fit$crf_control)
	    set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
      
	    log_obs_probs_by_state <- lapply(data_split, function(obs_dat) {
	      matrix(0, nrow = nrow(obs_dat$X), ncol = S)
	    })
	    
	    for(m in seq_len(cv_M[fit_ind])) {
	      temp <- lccrf_update_log_lik_parametric(betas = fit$component_fits[[fit_ind]]$beta_hat[[m]] / num_component_fits, omegas = omegas,
	                                                 active_cols = which(fit$component_fits[[fit_ind]]$active_cols[[m]]), beta_penalty_factor = 0, CRF_Xptr = CRF_Xptr,
	                                                 data_split = data_split, log_alpha = fit$component_fits[[fit_ind]]$log_alpha_hat[[m]],
	                                                 prev_log_obs_probs_by_state = log_obs_probs_by_state, crf_control = (fit$crf_control),
	                                                 retall = 1L)
	      log_obs_probs_by_state <- temp[[2]]
	    }
      
	    combined_log_class_probs <- CRF_log_lik_given_obs_probs(CRF_Xptr, data_split, log_obs_probs_by_state, time_scale = "marginal", retall = 1)[[2]]
      
	    for(dat_ind in seq_along(combined_log_class_probs)) {
	      y_pred_by_component_ind[[dat_ind]][, fit_ind] <- apply(combined_log_class_probs[[dat_ind]], 1, which.max)
	    }
	  }

    majority_vote_predictions <- lapply(y_pred_by_component_ind, function(ypci1) {
      as.integer(apply(ypci1, 1, function(ypc1) names(which.max(table(ypc1)))))
    })
    
    return(list(combined_predictions = majority_vote_predictions, component_predictions = y_pred_by_component_ind))
	} else {
		stop("Invalid prediction methods")
	}

#	if(identical(predict_method, "marginal")) {
#		predict_lccrf_given_obs_probs_marginal(CRF_Xptr, data_split, log_obs_probs_by_state)
#	} else if(identical(predict_method, "Viterbi")) {
#		predict_lccrf_given_obs_probs_Viterbi(CRF_Xptr, data_split, log_obs_probs_by_state)
#	} else {
#		stop("Invalid prediction method: must be marginal or Viterbi")
#	}
}

predict_lccrf_given_obs_probs_marginal <- function(CRF_Xptr, observedHMMs, log_obs_probs_by_state, return_class_probs = FALSE) {
	temp <- CRF_log_lik_given_obs_probs(CRF_Xptr, observedHMMs, log_obs_probs_by_state, time_scale = "marginal", retall = 1)
	if(return_class_probs) {
		return(list(predictions = lapply(temp[[2]], function(cp1) { apply(cp1, 1, which.max) }),
			log_class_probs = temp[[2]])
		)
	} else {
		return(lapply(temp[[2]], function(cp1) { apply(cp1, 1, which.max) }))
	}
}

predict_lccrf_given_obs_probs_Viterbi <- function(CRF_Xptr, observedHMMs, log_obs_probs_by_state) {
	stop("this method is not yet implemented")
}

calc_marginal_class_probs_given_log_obs_probs_CRF <- function(CRF_Xptr, observedHMMs, log_obs_probs_by_state, log) {
	temp <- CRF_log_lik_given_obs_probs(CRF_Xptr, observedHMMs, log_obs_probs_by_state, time_scale = "marginal", retall = 1)
	combined_log_ps <- rbind.fill.matrix(temp[[2]])
	if(log) {
		return(combined_log_ps)
	} else {
		return(exp(combined_log_ps))
	}
}



# C interface -- transition matrix should be set already
globalSeqCRF_update_log_lik_C_interface <- function(HMM_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
    betas, cols, retall) {
    .Call("parametric_CRF_update_log_lik_C", HMM_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        as.numeric(betas), as.integer(cols), as.integer(retall))
}

# interface with betas as first parameter -- sets transition matrix and calls C interface
globalSeqCRF_update_log_lik <- function(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, retall) {
    set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
    
    globalSeqCRF_update_log_lik_C_interface(CRF_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        betas, cols, retall)
}

globalSeqCRF_update_neg_log_lik <- function(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state) {
    return(-1 * globalSeqCRF_update_log_lik(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
            prev_log_obs_probs_by_state, retall = FALSE))
}

# interface with combined omegas and betas as first parameter -- sets transition matrix and calls C interface
globalSeqCRF_update_log_lik_combinedparamsfirst <- function(params, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, retall) {
    S <- ncol(prev_log_obs_probs_by_state[[1]])
    
    omega_hat <- params[seq_len(S^2 - 1)]
    omegas <- lapply(seq_len(S - 1), function(s) {
            omega_hat[(s - 1) * (S) + seq_len(S)]
        })
    omegas[[S]] <- omega_hat[(S - 1) * (S) + seq_len(S - 1)]
    set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
    
    betas <- params[seq(from = S^2, to = length(params), by = 1)]
    
    globalSeqCRF_update_log_lik_C_interface(CRF_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        betas, cols, retall)
}

globalSeqCRF_update_neg_log_lik_combinedparamsfirst <- function(params, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state) {
    return(-1 * globalSeqCRF_update_log_lik_combinedparamsfirst(params, cols, CRF_Xptr, obs_HMMs_split,
            prev_log_obs_probs_by_state, retall = 0L))
}





# C interface -- transition matrix should be set already
globalSeqCRF_update_gradient_log_lik_C_interface <- function(HMM_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
    betas, cols, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization) {
    .Call("parametric_CRF_update_gradient_log_lik_C", HMM_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        as.numeric(betas), as.integer(cols), as.integer(beta_diff_cols), as.integer(diff_wrt_omegas), as.integer(reduced_log_mat_parameterization))
}

# interface with betas as first parameter -- sets transition matrix and calls C interface
globalSeqCRF_update_gradient_log_lik <- function(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization) {
    set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
    
    globalSeqCRF_update_gradient_log_lik_C_interface(CRF_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        betas, cols, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization)
}

globalSeqCRF_update_gradient_neg_log_lik <- function(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization) {
    return(-1 * globalSeqCRF_update_gradient_log_lik(betas, omegas, cols, CRF_Xptr, obs_HMMs_split,
            prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization))
}

# interface with combined omegas and betas as first parameter -- sets transition matrix and calls C interface
globalSeqCRF_update_gradient_log_lik_combinedparamsfirst <- function(params, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization) {
    S <- ncol(prev_log_obs_probs_by_state[[1]])
    
    omega_hat <- params[seq_len(S^2 - 1)]
    omegas <- lapply(seq_len(S - 1), function(s) {
            omega_hat[(s - 1) * (S) + seq_len(S)]
        })
    omegas[[S]] <- omega_hat[(S - 1) * (S) + seq_len(S - 1)]
    set_trans_mat_from_omegas_CRF(CRF_Xptr, omegas)
    
    betas <- params[seq(from = S^2, to = length(params), by = 1)]
    
    globalSeqCRF_update_gradient_log_lik_C_interface(CRF_Xptr, obs_HMMs_split, prev_log_obs_probs_by_state,
        betas, cols, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization)
}

globalSeqCRF_update_gradient_neg_log_lik_combinedparamsfirst <- function(params, cols, CRF_Xptr, obs_HMMs_split,
    prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization) {
    return(-1 * globalSeqCRF_update_gradient_log_lik_combinedparamsfirst(params, cols, CRF_Xptr, obs_HMMs_split,
            prev_log_obs_probs_by_state, beta_diff_cols, diff_wrt_omegas, reduced_log_mat_parameterization))
}
