get_rot_mat_from_theta <- function(theta, D = 3) {
	cos_theta <- cos(theta)
	sin_theta <- sin(theta)

	if(D == 3) {
		rot_mat_1 <- diag(3)
		rot_mat_1[2, 2] <- cos_theta
		rot_mat_1[3, 3] <- cos_theta
		rot_mat_1[2, 3] <- -sin_theta
		rot_mat_1[3, 2] <- sin_theta
		rot_mat_2 <- diag(3)
		rot_mat_2[1, 1] <- cos_theta
		rot_mat_2[3, 3] <- cos_theta
		rot_mat_2[1, 3] <- sin_theta
		rot_mat_2[3, 1] <- -sin_theta
		rot_mat_3 <- diag(3)
		rot_mat_3[1, 1] <- cos_theta
		rot_mat_3[2, 2] <- cos_theta
		rot_mat_3[1, 2] <- -sin_theta
		rot_mat_3[2, 1] <- sin_theta

		rot_mat <- rot_mat_1 %*% rot_mat_2 %*% rot_mat_3
	} else if(D == 2) {
		rot_mat <- diag(2)
		rot_mat[1, 1] <- cos_theta
		rot_mat[2, 2] <- cos_theta
		rot_mat[1, 2] <- -sin_theta
		rot_mat[2, 1] <- sin_theta
	}

	return(rot_mat)
}

get_obs_dist_params <- function(obs_dist_normal, redundant_features_informative, bayes_error_rate_high, unequal_class_probs = FALSE) {
	if(obs_dist_normal) {
		if(! unequal_class_probs) {
			if(bayes_error_rate_high) {
#				group1_rot_angle <- list(0, pi/32, -pi/32)
				group1_rot_angle <- list(0, 0, 0)
				group1_offset <- list(0, 1.75, 2.0)
				group1_var_factor <- list(1, 0.9, 0.8)
	
#				group2_rot_angle <- list(0, pi/64, -pi/64)
				group2_rot_angle <- list(0, 0, 0)
				group2_offset <- list(rep(0, 3), rep(0.1, 3), rep(-0.1, 3))
#				group2_var_factor <- list(1, 0.9, 1.1)
				group2_var_factor <- list(1, 0.9, 0.9)
	
				group3_rot_angle <- list(0, 0, 0)
				group3_offset <- list(0, 0.6, -0.6)
				group3_var_factor <- list(1, 1, 1)
	
				group4_var_factor <- list(1, 1.2, 0.8)
			} else {
				group1_rot_angle <- list(0, pi/8, -pi/8)
				group1_offset <- list(0, 3.5, 3.5)
				group1_var_factor <- list(1, 0.6, 0.6)
	
#				group2_rot_angle <- list(pi/16, 0, -pi/16)
				group2_rot_angle <- list(pi/32, 0, -pi/32)
				group2_offset <- list(rep(0.25, 3), rep(0, 3), rep(-0.25, 3))
				group2_var_factor <- list(0.7, 1, 0.7)
	
				group3_rot_angle <- list(0, 0, 0)
				group3_offset <- list(0, 2.5, -2.5)
				group3_var_factor <- list(0.8, 1, 1)
	
				group4_var_factor <- list(1, 2, 0.5)
			}
		} else {
			if(bayes_error_rate_high) {
#				group1_rot_angle <- list(0, pi/32, -pi/32)
				group1_rot_angle <- list(0, 0, 0)
				group1_offset <- list(0, 0.2, 2.0)
				group1_var_factor <- list(1, 0.9, 0.8)
	
#				group2_rot_angle <- list(0, pi/64, -pi/64)
				group2_rot_angle <- list(0, 0, 0)
				group2_offset <- list(rep(0, 3), rep(0.05, 3), rep(-0.1, 3))
#				group2_var_factor <- list(1, 0.9, 1.1)
				group2_var_factor <- list(1, 1, 0.9)
	
				group3_rot_angle <- list(0, 0, 0)
				group3_offset <- list(0, 0.1, -0.6)
				group3_var_factor <- list(1, 1.1, 1)
	
				group4_var_factor <- list(2, 1.2, 0.8)
			} else {
				group1_rot_angle <- list(0, pi/8, -pi/8)
				group1_offset <- list(0, 3.5, 3.5)
				group1_var_factor <- list(1, 0.6, 0.6)
	
#				group2_rot_angle <- list(pi/16, 0, -pi/16)
				group2_rot_angle <- list(pi/32, 0, -pi/32)
				group2_offset <- list(rep(0.25, 3), rep(0, 3), rep(-0.25, 3))
				group2_var_factor <- list(0.7, 1, 0.7)
	
				group3_rot_angle <- list(0, 0, 0)
				group3_offset <- list(0, 2.5, -2.5)
				group3_var_factor <- list(0.8, 1, 1)
	
				group4_var_factor <- list(1, 2, 0.5)
			}
		}

		obs_dist_params_train <- lapply(seq_len(3), function(s) {
			if(s == 1) {
				### group 1 -- multimodality, group dimension = 3
				g1s1var <- 3
				g1s1covar <- 2.3
				group1 <- list(rho = 1,
					mus = list(c(0, 0, 0)),
					Sigmas = list(matrix(c(g1s1var, g1s1covar, g1s1covar, g1s1covar, g1s1var, g1s1covar, g1s1covar, g1s1covar, g1s1var), nrow = 3, ncol = 3)))

				### group 2 -- multimodality, overall non-linear shape, group dimension = 3
				g2s1m1var1 <- 4
				g2s1m1var2 <- 7
				g2s1m1var3 <- 4
				g2s1m1covar12 <- 3
				g2s1m1covar13 <- 3
				g2s1m1covar23 <- 0

				g2s1m2var1 <- 7
				g2s1m2var2 <- 7
				g2s1m2var3 <- 4
				g2s1m2covar12 <- -5.5
				g2s1m2covar13 <- -3
				g2s1m2covar23 <- 0

				g2s1m3var1 <- 9
				g2s1m3var2 <- 3
				g2s1m3var3 <- 4
				g2s1m3covar12 <- 3.5
				g2s1m3covar13 <- -3
				g2s1m3covar23 <- 0

				rot_mat <- get_rot_mat_from_theta(group2_rot_angle[[s]])

				group2 <- list(rho = c(0.3, 0.4, 0.3),
					mus = list(as.vector(rot_mat %*% matrix(c(-9, -1, 0))) + group2_offset[[s]], rep(0, 3) + group2_offset[[s]], as.vector(rot_mat %*% matrix(c(-1, -7, 0))) + group2_offset[[s]]),
					Sigmas = list(group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s1m1var1, g2s1m1covar12, g2s1m1covar13, g2s1m1covar12, g2s1m1var2, g2s1m1covar23, g2s1m1covar13, g2s1m1covar23, g2s1m1var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s1m2var1, g2s1m2covar12, g2s1m2covar13, g2s1m2covar12, g2s1m2var2, g2s1m2covar23, g2s1m2covar13, g2s1m2covar23, g2s1m2var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s1m3var1, g2s1m3covar12, g2s1m3covar13, g2s1m3covar12, g2s1m3var2, g2s1m3covar23, g2s1m3covar13, g2s1m3covar23, g2s1m3var3), nrow = 3, ncol = 3) %*% t(rot_mat)))

				### group 3 -- truncated normals, group dimension = 2
				g3s1var1 <- 8
				g3s1var2 <- 8
				g3s1covar <- 3

				rot_mat <- get_rot_mat_from_theta(group3_rot_angle[[s]], D = 2)

				group3 <- list(rho = 1,
					mus = list(as.vector(rot_mat %*% matrix(rep(group3_offset[[s]], 2)))),
					Sigmas = list(rot_mat %*% matrix(c(g3s1var1, g3s1covar, g3s1covar, g3s1var2), nrow = 2, ncol = 2) %*% t(rot_mat)))

				### group 4 -- same mean, different covariance
				group4 <- list(rho = 1,
					mus = list(c(0, 0)),
					Sigmas = list(group4_var_factor[[s]] * matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)))

			} else if(s == 2) {
				### group 1 -- multimodality, group dimension = 3
				g1s2var <- 1
				g1s2covar <- 0.5
				rot_mat <- get_rot_mat_from_theta(group1_rot_angle[[s]])
					group1 <- list(rho = c(0.25, 0.5, 0.25),
					mus = list(as.vector(rot_mat %*% matrix(rep(-group1_offset[[s]], 3))), rep(0, 3), as.vector(rot_mat %*% matrix(rep(group1_offset[[s]], 3)))),
					Sigmas = list(group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var), nrow = 3, ncol = 3) %*% t(rot_mat),
						group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var), nrow = 3, ncol = 3) %*% t(rot_mat),
						group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var, g1s2covar, g1s2covar, g1s2covar, g1s2var), nrow = 3, ncol = 3) %*% t(rot_mat)))
					
				### group 3 -- multimodality, overall non-linear shape, group dimension = 3
				g2s2m1var1 <- 4
				g2s2m1var2 <- 7
				g2s2m1var3 <- 4
				g2s2m1covar12 <- 3
				g2s2m1covar13 <- 3
				g2s2m1covar23 <- 0

				g2s2m2var1 <- 7
				g2s2m2var2 <- 7
				g2s2m2var3 <- 4
				g2s2m2covar12 <- -5.5
				g2s2m2covar13 <- -3
				g2s2m2covar23 <- 0

				g2s2m3var1 <- 9
				g2s2m3var2 <- 3
				g2s2m3var3 <- 4
				g2s2m3covar12 <- 3.5
				g2s2m3covar13 <- -3
				g2s2m3covar23 <- 0

				rot_mat <- get_rot_mat_from_theta(group2_rot_angle[[s]])

				group2 <- list(rho = c(0.3, 0.4, 0.3),
					mus = list(as.vector(rot_mat %*% matrix(c(-9, -1, 0))) + group2_offset[[s]], rep(0, 3) + group2_offset[[s]], as.vector(rot_mat %*% matrix(c(-1, -7, 0))) + group2_offset[[s]]),
					Sigmas = list(group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s2m1var1, g2s2m1covar12, g2s2m1covar13, g2s2m1covar12, g2s2m1var2, g2s2m1covar23, g2s2m1covar13, g2s2m1covar23, g2s2m1var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s2m2var1, g2s2m2covar12, g2s2m2covar13, g2s2m2covar12, g2s2m2var2, g2s2m2covar23, g2s2m2covar13, g2s2m2covar23, g2s2m2var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s2m3var1, g2s2m3covar12, g2s2m3covar13, g2s2m3covar12, g2s2m3var2, g2s2m3covar23, g2s2m3covar13, g2s2m3covar23, g2s2m3var3), nrow = 3, ncol = 3) %*% t(rot_mat)))

				### group 3 -- truncated normals, group dimension = 2
				g3s2var1 <- 8
				g3s2var2 <- 8
				g3s2covar <- 3

				rot_mat <- get_rot_mat_from_theta(group3_rot_angle[[s]], D = 2)

				group3 <- list(rho = 1,
					mus = list(as.vector(rot_mat %*% matrix(rep(group3_offset[[s]], 2)))),
					Sigmas = list(rot_mat %*% matrix(c(g3s2var1, g3s2covar, g3s2covar, g3s2var2), nrow = 2, ncol = 2) %*% t(rot_mat)))

				### group 4 -- same mean, different covariance
				group4 <- list(rho = 1,
					mus = list(c(0, 0)),
					Sigmas = list(group4_var_factor[[s]] * matrix(c(1, 0.4, 0.4, 1), nrow = 2, ncol = 2)))
			} else if(s == 3) {
				### group 1 -- multimodality, group dimension = 3
				g1s3var <- 1
				g1s3covar <- 0.5
				rot_mat <- get_rot_mat_from_theta(group1_rot_angle[[s]])

				group1 <- list(rho = c(0.25, 0.5, 0.25),
					mus = list(as.vector(rot_mat %*% matrix(rep(-group1_offset[[s]], 3))), rep(0, 3), as.vector(rot_mat %*% matrix(rep(group1_offset[[s]], 3)))),
					Sigmas = list(group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var), nrow = 3, ncol = 3) %*% t(rot_mat),
						group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var), nrow = 3, ncol = 3) %*% t(rot_mat),
						group1_var_factor[[s]] * rot_mat %*% matrix(c(g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var, g1s3covar, g1s3covar, g1s3covar, g1s3var), nrow = 3, ncol = 3) %*% t(rot_mat)))


				### group 2 -- multimodality, overall non-linear shape, group dimension = 3
				g2s3m1var1 <- 4
				g2s3m1var2 <- 7
				g2s3m1var3 <- 4
				g2s3m1covar12 <- 3
				g2s3m1covar13 <- 3
				g2s3m1covar23 <- 0

				g2s3m2var1 <- 7
				g2s3m2var2 <- 7
				g2s3m2var3 <- 4
				g2s3m2covar12 <- -5.5
				g2s3m2covar13 <- -3
				g2s3m2covar23 <- 0

				g2s3m3var1 <- 9
				g2s3m3var2 <- 3
				g2s3m3var3 <- 4
				g2s3m3covar12 <- 3.5
				g2s3m3covar13 <- -3
				g2s3m3covar23 <- 0

				rot_mat <- get_rot_mat_from_theta(group2_rot_angle[[s]])

				group2 <- list(rho = c(0.3, 0.4, 0.3),
					mus = list(as.vector(rot_mat %*% matrix(c(-9, -1, 0))) + group2_offset[[s]], rep(0, 3) + group2_offset[[s]], as.vector(rot_mat %*% matrix(c(-1, -7, 0))) + group2_offset[[s]]),
					Sigmas = list(group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s3m1var1, g2s3m1covar12, g2s3m1covar13, g2s3m1covar12, g2s3m1var2, g2s3m1covar23, g2s3m1covar13, g2s3m1covar23, g2s3m1var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s3m2var1, g2s3m2covar12, g2s3m2covar13, g2s3m2covar12, g2s3m2var2, g2s3m2covar23, g2s3m2covar13, g2s3m2covar23, g2s3m2var3), nrow = 3, ncol = 3) %*% t(rot_mat),
						group2_var_factor[[s]] * rot_mat %*% matrix(c(g2s3m3var1, g2s3m3covar12, g2s3m3covar13, g2s3m3covar12, g2s3m3var2, g2s3m3covar23, g2s3m3covar13, g2s3m3covar23, g2s3m3var3), nrow = 3, ncol = 3) %*% t(rot_mat)))

				### group 3 -- truncated normals, group dimension = 2
				g3s3var1 <- 8
				g3s3var2 <- 8
				g3s3covar <- 3

				rot_mat <- get_rot_mat_from_theta(group3_rot_angle[[s]], D = 2)

				group3 <- list(rho = 1,
					mus = list(as.vector(rot_mat %*% matrix(rep(group3_offset[[s]], 2)))),
					Sigmas = list(rot_mat %*% matrix(c(g3s3var1, g3s3covar, g3s3covar, g3s3var2), nrow = 2, ncol = 2) %*% t(rot_mat)))

				### group 4 -- same mean, different covariance
				group4 <- list(rho = 1,
					mus = list(c(0, 0)),
					Sigmas = list(group4_var_factor[[s]] * matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)))
			} else {
				stop("Invalid value of s")
			}

			obs_dist_params_one_state <- list(feature_group_params = list(group1 = group1, group2 = group2, group3 = group3, group4 = group4),
				obs_dist_normal = obs_dist_normal)

			if(redundant_features_informative) {
				obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
					diag(20)[, 1:10])
			} else {
				obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
			}

			return(obs_dist_params_one_state)
		})

		return(obs_dist_params_train)
	} else { # obs dist non normal
		if(bayes_error_rate_high) {
			obs_dist_params <- lapply(seq_len(3), function(s) {
				if(s == 1) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4 * c(0.3, 0.4, 0.3), 0.1, 0.4),
#								rho = 1,
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(0.9, 0.001),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(101),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0.001),
											gamma_beta = c(2, 0, 0.001)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(102),
											gamma_alpha = c(1.2),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0.001, 0),
											gamma_beta = c(2, 0, 0.002)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(98),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(98, 0),
											gamma_alpha = c(5, -0.01),
											gamma_beta = c(2, 0)
										),
										var3 = list(
											gamma_loc = c(98, 0, 0),
											gamma_alpha = c(4, 0.001, 0.001),
											gamma_beta = c(2, 0.001, 0.001)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(4),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(2.2, 0.001),
											gamma_beta = c(1, 0.001)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(2.1, 0.001, 0),
											gamma_beta = c(1.8, 0, 0)
										)
									)
								)
							),
							group2 = list(
								rho = rep(1/5, 5),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0.001, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(2.1)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0.001)
										),
										var3 = list(
											gamma_loc = c(160, -0.5, 0),
											gamma_alpha = c(1, 0.004, 0),
											gamma_beta = c(1, 0.001, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(7, 0.001),
											gamma_beta = c(0.52, 0)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.0, 0),
											gamma_beta = c(0.5, 0.00, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(6, 0.01),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.002, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(7.1, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									)
								)
							),
							group3 = list(
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(4.9),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-98.9, 1),
											gamma_beta = c(1, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(8),
											gamma_beta = c(1.1)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(4, 0.15),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(4, 0.002),
											gamma_beta = c(0.85, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(4.1),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(3.8, 0.001),
											gamma_beta = c(0.8, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = c(0.3, 0.3, 0.4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(10),
											gamma_beta = c(0.5)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(2, 0),
											gamma_beta = c(5, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(10),
											gamma_beta = c(5)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(9),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(40, 0.5),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1, 0)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)

					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}
						return(obs_dist_params_one_state)
				} else if(s == 2) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4 * c(0.3, 0.4, 0.3), 0.1, 0.4),
#								rho = 1,
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0.001)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0.002),
											gamma_beta = c(2, 0, 0.002)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(101),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0.002),
											gamma_beta = c(2, 0, 0.002)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(102),
											gamma_alpha = c(1.1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0.002, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(97),
											gamma_alpha = c(4.5),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(97, 0),
											gamma_alpha = c(4, 0.002),
											gamma_beta = c(2, 0)
										),
										var3 = list(
											gamma_loc = c(97, 0, 0),
											gamma_alpha = c(4, 0, 0),
											gamma_beta = c(2, 0.001, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(4),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(2, 0.001),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(2, 0, 0),
											gamma_beta = c(2, 0.002, 0)
										)
									)
								)
							),
							group2 = list(
#								rho = c(0, 1),
#								rho = c(0.5, 0.5),
								rho = rep(1/5, 5),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0.001)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.00),
											gamma_beta = c(1, 0.00)
										),
										var3 = list(
											gamma_loc = c(160, -0.5, 0),
											gamma_alpha = c(1, 0.002, 0),
											gamma_beta = c(1, 0.00, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									)
								)
							),
							group3 = list(
#								rho = c(0.5, 0.5),
#								rho = c(1),
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-99, 1),
											gamma_beta = c(1, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(8),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(4, 0.1),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(4, 0),
											gamma_beta = c(0.8, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(4, 0),
											gamma_beta = c(0.8, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = rep(1/3, 3),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(9.5),
											gamma_beta = c(0.5)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(2, 0.001),
											gamma_beta = c(5, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(11),
											gamma_beta = c(5)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1.1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(9),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(40, 0.5),
											gamma_alpha = c(8, 0.001),
											gamma_beta = c(1, 0.001)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)
					
					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}

					return(obs_dist_params_one_state)
				} else if(s == 3) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4 * c(0.3, 0.4, 0.3), 0.1, 0.4),
#								rho = 1,
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(1.5),
											gamma_beta = c(0.7)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(101),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(102),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(97),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(97, 0),
											gamma_alpha = c(4, -0.001),
											gamma_beta = c(2, 0)
										),
										var3 = list(
											gamma_loc = c(97, 0, 0),
											gamma_alpha = c(4, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(4),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(2, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(2, 0.001, 0.001),
											gamma_beta = c(2, 0, 0)
										)
									)
								)
							),
							group2 = list(
#								rho = c(0, 1),
#								rho = c(0.5, 0.5),
								rho = rep(1/5, 5),
#								rho = c(0.15, 0.25, 0.2, 0.15, 0.25),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0.000)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.00),
											gamma_beta = c(1, 0.00)
										),
										var3 = list(
											gamma_loc = c(160, -0.5, 0),
											gamma_alpha = c(1, 0.002, 0),
											gamma_beta = c(1, 0.00, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(7, 0.002),
											gamma_beta = c(0.5, 0.001)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									)
								)
							),
							group3 = list(
#								rho = c(0.5, 0.5),
#								rho = c(1),
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(0.8)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-99, 1),
											gamma_beta = c(1, 0.001)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(9),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(3.5, 0.1),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(4, 0),
											gamma_beta = c(1, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(104.5),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(4, 0),
											gamma_beta = c(0.8, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = rep(1/3, 3),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(10.5),
											gamma_beta = c(0.5)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(2.1, 0),
											gamma_beta = c(5, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(10),
											gamma_beta = c(5)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(8, 0.001),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(9),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(40, 0.5),
											gamma_alpha = c(8, 0.001),
											gamma_beta = c(1, 0)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)

					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}

					return(obs_dist_params_one_state)
				}
			})

#			normal_obs_dist_params <- get_obs_dist_params(obs_dist_normal = TRUE, redundant_features_informative = redundant_features_informative, bayes_error_rate_high = bayes_error_rate_high)
#			obs_dist_params <- convert_normal_obs_dist_params_to_nonnormal(normal_obs_dist_params)
		} else {
			obs_dist_params <- lapply(seq_len(3), function(s) {
				if(s == 1) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4, 0.1, 0.4),
#								rho = c(0, 1, 0),
#								rho = c(0, 0, 1),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100.5),
											gamma_alpha = c(1),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(2, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(98),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(98, 0),
											gamma_alpha = c(4, 0.001),
											gamma_beta = c(2, 0.001)
										),
										var3 = list(
											gamma_loc = c(98, 0, 0),
											gamma_alpha = c(3.5, 0, 0.002),
											gamma_beta = c(1.8, 0.001, 0.001)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(5),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(3, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(3, 0.001, 0),
											gamma_beta = c(2, 0, 0)
										)
									)
								)
							),
							group2 = list(
								rho = rep(1/5, 5),
#								rho = c(0, 1),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2.5)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0.001, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(4),
											gamma_beta = c(3)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0.01)
										),
										var3 = list(
											gamma_loc = c(160, -0.5, 0),
											gamma_alpha = c(1, 0.01, 0),
											gamma_beta = c(1, 0.001, 0.001)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(5, 0),
											gamma_beta = c(0.5, 0.001)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.002, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(1.1)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(7, 0.01),
											gamma_beta = c(0.5, 0.002)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.02, 0),
											gamma_beta = c(0.5, 0.001, 0.001)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(9, -0.001),
											gamma_beta = c(0.5, 0.001)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(3, 0.02, 0),
											gamma_beta = c(0.7, 0.002, 0)
										)
									)
								)
							),
							group3 = list(
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(4.5),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-98.9, 1),
											gamma_beta = c(1, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(8),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(4.4, 0.01),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(4, 0),
											gamma_beta = c(0.9, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(3, 0),
											gamma_beta = c(0.8, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = c(0.3, 0.3, 0.4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(10),
											gamma_beta = c(0.5)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(2, 0),
											gamma_beta = c(5, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(10),
											gamma_beta = c(5)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(9),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(40, 0.5),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1, 0)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)

					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}
						return(obs_dist_params_one_state)
				} else if(s == 2) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4 * c(0.3, 0.4, 0.3), 0.1, 0.4),
#								rho = 1,
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(101),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(102),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(98),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(98, 0),
											gamma_alpha = c(4, 0),
											gamma_beta = c(2, 0)
										),
										var3 = list(
											gamma_loc = c(98, 0, 0),
											gamma_alpha = c(4, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(4),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(2, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(2, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									)
								)
							),
							group2 = list(
#								rho = c(0, 1),
#								rho = c(0.5, 0.5),
								rho = rep(1/5, 5),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0.001)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.00),
											gamma_beta = c(1, 0.00)
										),
										var3 = list(
											gamma_loc = c(160, -0.5, 0),
											gamma_alpha = c(1, 0.002, 0),
											gamma_beta = c(1, 0.00, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									)
								)
							),
							group3 = list(
#								rho = c(0.5, 0.5),
#								rho = c(1),
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5.5),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-99, 1),
											gamma_beta = c(1, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(1.5)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(4.5, 0.1),
											gamma_beta = c(1, 0.1)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2.5),
											gamma_beta = c(2.5)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(3.5, 0),
											gamma_beta = c(0.8, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(3, 0),
											gamma_beta = c(0.6, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = rep(1/3, 3),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(9),
											gamma_beta = c(0.6)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(2, 0.001),
											gamma_beta = c(5, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(11),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1.5, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(9),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(43, 0.48),
											gamma_alpha = c(8, 0),
											gamma_beta = c(1, 0.001)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)
					
					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}

					return(obs_dist_params_one_state)
				} else if(s == 3) {
					obs_dist_params_one_state <- list(
						feature_group_params = list(
							group1 = list(
								rho = c(0.4 * c(0.3, 0.4, 0.3), 0.1, 0.4),
#								rho = 1,
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(0.5)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1.5, 0.005),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2.5, 0.002, 0.002)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(101),
											gamma_alpha = c(1.1),
											gamma_beta = c(1.5)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(2, 0.01),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0.01),
											gamma_beta = c(2, 0, 0.01)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(102),
											gamma_alpha = c(1),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(1.5, 0),
											gamma_beta = c(1.5, 0)
										),
										var3 = list(
											gamma_loc = c(0, 1, 0),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(1, .01, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(97),
											gamma_alpha = c(5),
											gamma_beta = c(1.5)
										),
										var2 = list(
											gamma_loc = c(97, 0),
											gamma_alpha = c(4, -0.001),
											gamma_beta = c(2, 0.001)
										),
										var3 = list(
											gamma_loc = c(97, 0, 0),
											gamma_alpha = c(4, 0.001, 0.001),
											gamma_beta = c(2, 0, 0.001)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(104),
											gamma_alpha = c(4),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(210, -1),
											gamma_alpha = c(2, 0),
											gamma_beta = c(1, 0)
										),
										var3 = list(
											gamma_loc = c(80, 0.2, 0),
											gamma_alpha = c(2, 0.005, 0.005),
											gamma_beta = c(2, 0, 0)
										)
									)
								)
							),
							group2 = list(
#								rho = c(0, 1),
#								rho = c(0.5, 0.5),
#								rho = rep(1/5, 5),
								rho = c(0.15, 0.25, 0.2, 0.15, 0.25),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(1, 0.001),
											gamma_beta = c(1, 0.000)
										),
										var3 = list(
											gamma_loc = c(0, 0, 1),
											gamma_alpha = c(1, 0, 0),
											gamma_beta = c(2, 0, 0)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0.1),
											gamma_alpha = c(1, 0.00),
											gamma_beta = c(1, 0.00)
										),
										var3 = list(
											gamma_loc = c(153, -0.45, 0),
											gamma_alpha = c(1, 0.002, 0),
											gamma_beta = c(1, 0.00, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(105),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(55, 0.4),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(215, -1, 0),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(115),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(105, 0),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(325, -1, -1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									),
									component5 = list(
										var1 = list(
											gamma_loc = c(125),
											gamma_alpha = c(15),
											gamma_beta = c(0.9)
										),
										var2 = list(
											gamma_loc = c(145, -0.3),
											gamma_alpha = c(7, 0),
											gamma_beta = c(0.5, 0)
										),
										var3 = list(
											gamma_loc = c(-25, 0, 1),
											gamma_alpha = c(5, 0.01, 0),
											gamma_beta = c(0.5, 0.001, 0)
										)
									)
								)
							),
							group3 = list(
#								rho = c(0.5, 0.5),
#								rho = c(1),
								rho = rep(1/4, 4),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(5),
											gamma_beta = c(0.8)
										),
										var2 = list(
											gamma_loc = c(0, 1),
											gamma_alpha = c(-99, 1),
											gamma_beta = c(1, 0.001)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(9),
											gamma_beta = c(1)
										),
										var2 = list(
											gamma_loc = c(190, -1),
											gamma_alpha = c(3.5, 0.1),
											gamma_beta = c(1, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(110),
											gamma_alpha = c(2),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(90, 0),
											gamma_alpha = c(4, 0),
											gamma_beta = c(1, 0)
										)
									),
									component4 = list(
										var1 = list(
											gamma_loc = c(104.5),
											gamma_alpha = c(4),
											gamma_beta = c(2)
										),
										var2 = list(
											gamma_loc = c(-35, 1),
											gamma_alpha = c(4, 0),
											gamma_beta = c(0.8, 0.01)
										)
									)
								)
							),
							group4 = list(
								rho = rep(1/3, 3),
								component_params = list(
									component1 = list(
										var1 = list(
											gamma_loc = c(100),
											gamma_alpha = c(7),
											gamma_beta = c(0.7)
										),
										var2 = list(
											gamma_loc = c(100, 0),
											gamma_alpha = c(3, 0.005),
											gamma_beta = c(4, 0.01)
										)
									),
									component2 = list(
										var1 = list(
											gamma_loc = c(80),
											gamma_alpha = c(8),
											gamma_beta = c(6)
										),
										var2 = list(
											gamma_loc = c(120, 0),
											gamma_alpha = c(7, 0.01),
											gamma_beta = c(1.5, 0)
										)
									),
									component3 = list(
										var1 = list(
											gamma_loc = c(90),
											gamma_alpha = c(10),
											gamma_beta = c(4)
										),
										var2 = list(
											gamma_loc = c(40, 0.5),
											gamma_alpha = c(5, 0.3),
											gamma_beta = c(1, 0.1)
										)
									)
								)
							)
						),
						obs_dist_normal = FALSE
					)

					if(redundant_features_informative) {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- cbind(diag(20) + cbind(diag(2, nrow = 20), c(2, rep(0, 19)))[, 2:21],
							diag(20)[, 1:10])
					} else {
						obs_dist_params_one_state$redundant_vars_lc_matrix <- matrix(0, nrow = 20, ncol = 30)
					}

					return(obs_dist_params_one_state)
				}
			})

#			normal_obs_dist_params <- get_obs_dist_params(obs_dist_normal = TRUE, redundant_features_informative = redundant_features_informative, bayes_error_rate_high = bayes_error_rate_high)
#			obs_dist_params <- convert_normal_obs_dist_params_to_nonnormal(normal_obs_dist_params)
		}

		return(obs_dist_params)
	}
}



robs_dist <- function(n, obs_dist_params_one_state) {
	if(n == 0) {
		return(numeric(0))
	}
	
	if(obs_dist_params_one_state$obs_dist_normal) {
		rMM <- function(n, obs_dist_params_one_state_and_var_group) {
			return( rGMM(n, rho = obs_dist_params_one_state_and_var_group$rho, mus = obs_dist_params_one_state_and_var_group$mus, Sigmas = obs_dist_params_one_state_and_var_group$Sigmas) )
		}
	} else {
		rMM <- function(n, obs_dist_params_one_state_and_var_group) {
			return( rGammaMM(n, rho = obs_dist_params_one_state_and_var_group$rho, component_params = obs_dist_params_one_state_and_var_group$component_params) )
		}
	}

	result <- cbind(rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group1),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group2),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group3),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group4),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group1),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group2),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group3),
		rMM(n, obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group4))

	redundant_component_means <- result %*% obs_dist_params_one_state$redundant_vars_lc_matrix
	result <- cbind(result, redundant_component_means + rmvnorm(n,
		mean = rep(0, ncol(obs_dist_params_one_state$redundant_vars_lc_matrix)),
		sigma = diag(ncol(obs_dist_params_one_state$redundant_vars_lc_matrix))
	))

	return(result)
}

dobs_dist <- function(X, obs_dist_params_one_state, log) {
	if(! is.matrix(X) ) {
		X <- matrix(X, nrow = 1)
	}

	if(obs_dist_params_one_state$obs_dist_normal) {
		dMM <- function(X, obs_dist_params_one_state_and_var_group, log) {
			return( dGMM(X, rho = obs_dist_params_one_state_and_var_group$rho, mus = obs_dist_params_one_state_and_var_group$mus, Sigmas = obs_dist_params_one_state_and_var_group$Sigmas, log = log) )
		}
	} else {
		dMM <- function(X, obs_dist_params_one_state_and_var_group, log) {
			return( dGammaMM(X, rho = obs_dist_params_one_state_and_var_group$rho, component_params = obs_dist_params_one_state_and_var_group$component_params, log = log) )
		}
	}

	log_result <- dMM(X[, 1:3, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group1, log = TRUE) +
		dMM(X[, 4:6, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group2, log = TRUE) +
		dMM(X[, 7:8, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group3, log = TRUE) +
		dMM(X[, 9:10, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group4, log = TRUE) +
		dMM(X[, 11:13, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group1, log = TRUE) +
		dMM(X[, 14:16, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group2, log = TRUE) +
		dMM(X[, 17:18, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group3, log = TRUE) +
		dMM(X[, 19:20, drop = FALSE], obs_dist_params_one_state_and_var_group = obs_dist_params_one_state$feature_group_params$group4, log = TRUE)

	redundant_component_means <- X[, 1:20, drop = FALSE] %*% obs_dist_params_one_state$redundant_vars_lc_matrix
	log_result <- log_result + dmvnorm(X[, 21:50, drop = FALSE] - redundant_component_means, mean = rep(0, 30), sigma = diag(30), log = TRUE)

	if(log) {
		return(log_result)
	} else {
		return(exp(log_result))
	}
}

#robs_dist <- function(n, obs_dist_params_one_state) {
#	if(! obs_dist_params_one_state$obs_dist_complex) {
#		return(robs_dist_GMM(n, obs_dist_params_one_state))
#	} else {
#		return(robs_dist_gamma(n, obs_dist_params_one_state))
#	}
#}


robs_dist_GMM <- function(n, obs_dist_params_one_state) {
	if(obs_dist_params_one_state$D == 10) {
		return(cbind(rGMM(n, obs_dist_params_one_state$component1$rho, obs_dist_params_one_state$component1$mus, obs_dist_params_one_state$component1$Sigmas),
			rTGMM(n, obs_dist_params_one_state$component2$rho, obs_dist_params_one_state$component2$mus, obs_dist_params_one_state$component2$Sigmas,
				obs_dist_params_one_state$component2$lower, obs_dist_params_one_state$component2$upper, obs_dist_params_one_state$component2$linearconst),
#			rGMM(n, obs_dist_params_one_state$component2$rho, obs_dist_params_one_state$component2$mus, obs_dist_params_one_state$component2$Sigmas),
			rGMM(n, obs_dist_params_one_state$component3$rho, obs_dist_params_one_state$component3$mus, obs_dist_params_one_state$component3$Sigmas),
			rGMM(n, obs_dist_params_one_state$component4$rho, obs_dist_params_one_state$component4$mus, obs_dist_params_one_state$component4$Sigmas)))
	} else if(obs_dist_params_one_state$D == 80) {
		X <- c()
		for(group_ind in seq(from = 0, to = 2)) {
			start_ind <- group_ind * 20 + 1
			X <- cbind(X, robs_dist_GMM(n, obs_dist_params_one_state[[group_ind + 1]]))

			for(comp_ind in seq(from = 1, length = 10)) {
				X <- cbind(X, rnorm(n, mean = X[, start_ind + comp_ind - 1], sd = obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[comp_ind]))
			}
		}

		X <- cbind(X, rmvnorm(n, mean = rep(0, 20), sigma = diag(20)))
	}
}

robs_dist_gamma <- function(n, obs_dist_params_one_state) {
	if(obs_dist_params_one_state$D == 10) {
		X <- c()
		for(group_params in obs_dist_params_one_state$feature_group_params) {
			X <- cbind(X, rGammaMM(n, group_params$rho, group_params$component_params))
		}
	} else if(obs_dist_params_one_state$D == 80) {
		X <- c()
		for(group_ind in seq_along(obs_dist_params_one_state$feature_group_params)) {
			start_ind <- group_ind * 20 + 1
			X <- cbind(X, robs_dist_gamma(n, obs_dist_params_one_state[[group_ind + 1]]))

			for(comp_ind in seq(from = 1, length = 10)) {
				X <- cbind(X, rnorm(n, mean = X[, start_ind + comp_ind - 1], sd = obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[comp_ind]))
			}
		}

		X <- cbind(X, rmvnorm(n, mean = rep(0, 20), sigma = diag(20)))
	}

	return(X)
}




dobs_dist_gamma <- function(X, obs_dist_params_one_state, log) {
	if(!is.matrix(X)) {
		X <- matrix(X, nrow = 1)
	}
	if(obs_dist_params_one_state$D == 10) {
		log_result <- 0
		cum_num_vars <- 0
		for(group_params in obs_dist_params_one_state$feature_group_params) {
			log_result <- log_result + dGammaMM(X[, cum_num_vars + seq_along(group_params$component_params$component1)], group_params$rho, group_params$component_params, log = TRUE)
			cum_num_vars <- cum_num_vars + length(group_params$component_params$component1)
		}
	} else if(obs_dist_params_one_state$D == 80) {
		log_result <- 0
		for(group_ind in seq(from = 0, to = 2)) {
			start_ind <- group_ind * 20 + 1
			log_result <- log_result + dobs_dist_gamma(X[seq(from = start_ind, length = 10)], obs_dist_params_one_state[[group_ind + 1]], log = TRUE)
			log_result <- log_result + dmvnorm(matrix(X[seq(from = start_ind + 10, length = 3)], nrow = 1), mean = X[seq(from = start_ind, length = 3)], sigma = diag(obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[seq(from = 1, length = 3)]^2), log = TRUE)
			log_result <- log_result + dmvnorm(matrix(X[seq(from = start_ind + 10 + 3, length = 2)], nrow = 1), mean = X[seq(from = start_ind + 3, length = 2)], sigma = diag(obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[seq(from = 4, length = 2)]^2), log = TRUE)
			log_result <- log_result + dmvnorm(matrix(X[seq(from = start_ind + 10 + 5, length = 3)], nrow = 1), mean = X[seq(from = start_ind + 5, length = 3)], sigma = diag(obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[seq(from = 6, length = 3)]^2), log = TRUE)
			log_result <- log_result + dmvnorm(matrix(X[seq(from = start_ind + 10 + 8, length = 2)], nrow = 1), mean = X[seq(from = start_ind + 8, length = 2)], sigma = diag(obs_dist_params_one_state[[group_ind + 1]]$comp_offset_sd[seq(from = 9, length = 2)]^2), log = TRUE)
		}
		log_result <- log_result + dmvnorm(matrix(X[seq(from = 61, to = 80)], nrow = 1), mean = rep(0, 20), sigma = diag(20), log = TRUE)
	}

	if(log) {
		return(log_result)
	} else {
		return(exp(log_result))
	}
}
