get_C_dbl_max <- function() {
	return(.Call("get_dbl_max"))
}

logspace_sub <- function(logx, logy) {
	return(.Call("logspace_sub_R_C_interface", as.numeric(logx), as.numeric(logy)))
}

logspace_add <- function(logx, logy) {
	return(.Call("logspace_add_R_C_interface", as.numeric(logx), as.numeric(logy)))
}

logspace_sum_matrix_rows <- function(logX) {
	return(.Call("logspace_sum_matrix_rows", as.numeric(logX), as.integer(nrow(logX)), as.integer(ncol(logX))))
}

logspace_sub_matrix_rows <- function(logX) {
	if(!is.matrix(logX) || !identical(ncol(logX), 2L))
		stop("logX must be a matrix with 2 columns")

	return(.Call("logspace_sub_matrix_rows", as.numeric(logX), as.integer(nrow(logX))))
}

dMVN <- function(X, mu, Sigma, log = FALSE) {
	stop("This code is unreliable")
	log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus
	cat(log_det_Sigma)
	cat("\n")
	
	.Call("dMVN_multiobs_R_interface", as.numeric(X), nrow(X), ncol(X), as.numeric(mu), as.numeric(Sigma), as.numeric(log_det_Sigma), as.integer(log))
}

dMVN_prec <- function(X, mu, prec, log = FALSE) {
	stop("This code is unreliable")
	log_det_Sigma <- -determinant(prec, logarithm = TRUE)$modulus
	cat(log_det_Sigma)
	cat("\n")
	
	.Call("dMVN_Prec_multiobs_R_interface", as.numeric(X), nrow(X), ncol(X), as.numeric(mu), as.numeric(prec), as.numeric(log_det_Sigma), as.integer(log))
}

dMVN_diagprec <- function(X, mu, prec, log = FALSE) {
	stop("This code is unreliable")
	log_det_Sigma <- -sum(log(prec))
	cat(log_det_Sigma)
	cat("\n")
	
	.Call("dMVN_DiagPrec_multiobs_R_interface", as.numeric(X), nrow(X), ncol(X), as.numeric(mu), as.numeric(prec), as.numeric(log_det_Sigma), as.integer(log))
}

dGMM <- function(x, rhos = rep(1/length(mus), length(mus)), mus, Prec = NULL, Sigmas = NULL, log_det_Sigmas = NULL, log = FALSE) {
	if(log) {
		retlog <- 1L
	} else {
		retlog <- 0L
	}
	
	if(!identical(length(rhos), length(mus)))
		stop("rhos and mus must have the same length")
	
	if(!is.null(Sigmas)) {
		if(is.matrix(Sigmas)) {
			if(is.null(log_det_Sigmas))
				log_det_Sigmas <- determinant(Sigmas, logarithm = TRUE)$modulus
			
			if(! is.matrix(x)) {
				warning("x is not a matrix -- converting to matrix now")
				x <- matrix(x, ncol = ncol(Sigmas))
			}

			res <- matrix(-Inf, nrow = nrow(x), ncol = length(rhos))
			for(component_ind in seq_along(rhos)) {
				res[, component_ind] <- log(rhos[component_ind]) + dmvnorm(x, mean = mus[[component_ind]], sigma = Sigmas[[component_ind]], log = TRUE)
			}
			res <- logspace_sum_matrix_rows(res)

			if(log) {
				return(res)
			} else {
				return(exp(res))
			}
			
#			return(.Call("dGMM_sameSigma_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
#							as.numeric(rhos), mus, as.numeric(Sigmas), as.numeric(log_det_Sigmas), length(mus), retlog))
		} else if(is.list(Sigmas)) {
			if(!identical(length(rhos), length(Sigmas)))	
				stop("If Sigmas is a list, it must have the same length as rhos and mus")
			
			if(is.null(log_det_Sigmas))
				log_det_Sigmas <- lapply(Sigmas, function(Sigma) determinant(Sigma, logarithm = TRUE)$modulus)
#			browser()
			
			if(! is.matrix(x)) {
				warning("x is not a matrix -- converting to matrix now")
				x <- matrix(x, ncol = ncol(Sigmas[[1]]))
			}

			res <- matrix(-Inf, nrow = nrow(x), ncol = length(rhos))
			for(component_ind in seq_along(rhos)) {
				res[, component_ind] <- log(rhos[component_ind]) + dmvnorm(x, mean = mus[[component_ind]], sigma = Sigmas[[component_ind]], log = TRUE)
			}
			res <- logspace_sum_matrix_rows(res)

			if(log) {
				return(res)
			} else {
				return(exp(res))
			}
			
#			return(.Call("dGMM_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
#							as.numeric(rhos), mus, Sigmas, log_det_Sigmas, length(mus), retlog))
		}
	} else if(!is.null(Prec)) {
		stop("This code is unreliable")
		if(is.vector(Prec) && is.numeric(Prec)) {
			if(!identical(length(Prec), length(mus[[1]]))) {
				stop("If Prec is a numeric vector, it must have the same length as each mean vector")
			}
			
			if(is.null(log_det_Sigmas)) {
				log_det_Sigmas <- -sum(log(Prec))
			}
			
			
			if(! is.matrix(x)) {
				warning("x is not a matrix -- converting to matrix now")
				x <- matrix(x, ncol = length(Prec))
			}
			
			return(.Call("dGMM_sameDiagPrec_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
							as.numeric(rhos), mus, as.numeric(Prec), as.numeric(log_det_Sigmas), length(mus), retlog))
		} else if(is.matrix(Prec)) {
			if(is.null(log_det_Sigmas))
				log_det_Sigmas <- -determinant(Prec, logarithm = TRUE)$modulus
			
			
			if(! is.matrix(x)) {
				warning("x is not a matrix -- converting to matrix now")
				x <- matrix(x, ncol = ncol(Prec))
			}
			
			return(.Call("dGMM_samePrec_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
							as.numeric(rhos), mus, as.numeric(Prec), as.numeric(log_det_Sigmas), length(mus), retlog))
		} else if(is.list(Prec)) {
			if(!identical(length(rhos), length(Prec)))	
				stop("If Prec is a list, it must have the same length as rhos and mus")
			
			if(is.vector(Prec[[1]])) {
				if(!identical(length(Prec[[1]]), length(mus[[1]]))) {
					stop("If Prec is a vector, it must have the same length as each mean vector")
				}
				
				if(is.null(log_det_Sigmas)) {
					log_det_Sigmas <- lapply(Prec, function(P) -sum(log(P)))
				}
				
				
				if(! is.matrix(x)) {
					warning("x is not a matrix -- converting to matrix now")
					x <- matrix(x, ncol = length(Prec[[1]]))
				}
				
				return(.Call("dGMM_DiagPrec_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
								as.numeric(rhos), mus, Prec, log_det_Sigmas, length(mus), retlog))
			} else if(is.matrix(Prec[[1]])) {
				if(is.null(log_det_Sigmas))
					log_det_Sigmas <- lapply(Prec, function(P) -determinant(P, logarithm = TRUE)$modulus)
				
				if(! is.matrix(x)) {
					warning("x is not a matrix -- converting to matrix now")
					x <- matrix(x, ncol = ncol(Prec[[1]]))
				}
				
				return(.Call("dGMM_Prec_R_Interface", as.numeric(x), as.integer(nrow(x)), as.integer(ncol(x)),
								as.numeric(rhos), mus, Prec, log_det_Sigmas, length(mus), retlog))
			} else {
				stop("If Prec is a list, it must be a list of vectors or matrices.")
			}
		}
	} else {
		stop("Sigmas or Prec must be either a list of Sigma/Prec matrices for each of the M components, or a single matrix to be used for all components.")
	}
}
