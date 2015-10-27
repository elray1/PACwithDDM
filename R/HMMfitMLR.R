fit_MLR <- function(y, X, offsets, weights = rep(1, nrow(X))) {
  n_obs_classes <- length(unique(y))
  
  if(identical(1L, n_obs_classes)) {
    stop("y contains only one class")
#    result <- rep(-10, S * num_active_cols)
#    class_ind <- (levels(obs_HMMs_concat$y) == unique(obs_HMMs_concat$y))
#    result[class_ind - 1 + seq_len(num_active_cols)] <- 0
  } else {
    if(!is.factor(y))
      y <- factor(y)
    
    S <- length(levels(y))
    
    reduced_levels <- levels(y)[levels(y) %in% y]
    reduced_level_inds <- which(levels(y) %in% reduced_levels)
    y_red <- factor(y, levels = reduced_levels)
    
    y_for_mlr <- matrix(0, length(y_red), n_obs_classes)
    y_for_mlr[cbind(seq_along(y_red), sapply(y_red, function(yy) which(reduced_levels == yy)))] <- 1
    
    if(isTRUE(all.equal(X[, 1], rep(1, nrow(X)))) && ncol(X) > 1)
      X <- X[, 2:ncol(X)]
    
    D <- ncol(X) + 1
    
    if(!missing(offsets)) {
      temp <- multinom(y_for_mlr ~ X + offset(offsets), weights = weights)
    } else {
      temp <- multinom(y_for_mlr ~ X, weights = weights)
    }
    
    result <- rep(-Inf, S * D)
    
    inds <- c()
    for(rli in reduced_level_inds)
      inds <- c(inds, (rli - 1) * D + seq_len(D))
    
    result[inds] <- c(rep(0, D), as.numeric(t(coef(temp))))
  }

  return(list(S = S, D = D, coef = result))
}

calc_classprobs_MLR <- function(X, MLR_fit, log = FALSE) {
  if(!isTRUE(all.equal(X[, 1], rep(1, nrow(X)))) && ncol(X) > 1)
    X <- cbind(rep(1, nrow(X)), X)
  
  if(ncol(X) != MLR_fit$D)
    stop("Incorrect number of columns in X")
  
  lcp <- sapply(seq_len(MLR_fit$S), function(s) {
    betas <- matrix(MLR_fit$coef[(s - 1) * MLR_fit$D + seq_len(MLR_fit$D)])
    X %*% betas
  })
  
  normconst <- logspace_sum_matrix_rows(lcp)
  lcp <- sweep(lcp, 1, normconst, `-`)
  
  if(log) {
    return(lcp)
  } else {
    return(exp(lcp))
  }
}

predict_MLR <- function(X, MLR_fit) {
  classprobs <- calc_classprobs_MLR(X, MLR_fit, log = TRUE)
  return(apply(classprobs, 1, which.max))
}
