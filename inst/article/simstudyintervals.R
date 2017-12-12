<<fitSimStudyResultsLMEModel, cache = FALSE, echo = FALSE>>=
suppressMessages(suppressWarnings(library("nlme")))
suppressMessages(suppressWarnings(library("multcomp")))
# fit mixed effects model with a separate mean for each combination of
# data set, accelerometer location, response variable and fit method, random effect for subject, and
# variance specific to combination of subject, data set, response variable, and location


# results_fit <- lme(fixed = prop_correct ~ location * fit_method * data_set * response,
# 	random = ~ 1 | subject,
# #	weights = varIdent(form = ~ 1 | subject * location),
# #	weights = varIdent(form = ~ 1 | subject * location * data_set * response),
# 	weights = varIdent(form = ~ 1 | location * fit_method * data_set * response),
# 	data = combined_results_by_subject,
# 	control = lmeControl(maxIter = 500, msMaxIter = 500, niterEM = 250, msMaxEval = 2000))

#table(paste(combined_results_by_subject$location, combined_results_by_subject$fit_method, combined_results_by_subject$data_set, combined_results_by_subject$response, combined_results_by_subject$subject))

results_fit <- lme(fixed = prop_correct ~ obs_dist_complex * time_dep * fit_method,
 	random = ~ 1 | sim_ind,
#	weights = varIdent(form = ~ 1 | subject * location),
#	weights = varIdent(form = ~ 1 | subject * location * data_set * response),
 	weights = varIdent(form = ~ 1 | obs_dist_complex * time_dep * fit_method),
 	data = simStudyResults,
 	control = lmeControl(maxIter = 500, msMaxIter = 500, niterEM = 250, msMaxEval = 2000))


# assemble a data frame with estimates of relevant linear combinations of parameters and CIs.
# We want estimates for:
#  - the mean for each combination of location and classification method
#  - the difference in means between each pair of methods within location
#  - the difference in means between each location within each method
unique_fit_methods <- as.character(unique(simStudyResults$fit_method))
unique_obs_dist_complex <- as.character(unique(simStudyResults$obs_dist_complex))
unique_time_dep <- as.character(unique(simStudyResults$time_dep))

num_fit_methods <- length(unique_fit_methods)
num_obs_dist_complex <- length(unique_obs_dist_complex)
num_time_dep <- length(unique_time_dep)

unique_fit_method_descriptors <- paste0("fit_method", sort(unique_fit_methods))
unique_obs_dist_complex_descriptors <- paste0("obs_dist_complex", sort(unique_obs_dist_complex))
unique_time_dep_descriptors <- paste0("time_dep", sort(unique_time_dep))

lc_df <- expand.grid(
  fit_method = unique_fit_methods,
  obs_dist_complex = unique_obs_dist_complex,
  time_dep = unique_time_dep,
  stringsAsFactors = FALSE)

lc_df$fit_method_descriptor <- paste0("fit_method", lc_df$fit_method)
lc_df$obs_dist_complex_descriptor <- paste0("obs_dist_complex", lc_df$obs_dist_complex)
lc_df$time_dep_descriptor <- paste0("time_dep", lc_df$time_dep)
lc_df$name <- apply(as.matrix(lc_df[, 1:4]), 1, paste, collapse = "-")

num_leading_cols <- ncol(lc_df)
coef_cols <- seq(
  from = num_leading_cols + 1,
  length = num_fit_methods * num_obs_dist_complex * num_time_dep
)

# corresponding indicator vector for each coefficient
coef_names <- names(fixef(results_fit))
unique_coef_name_component_descriptors <- unique(unlist(strsplit(coef_names, ":")))
intercept_fit_method <- unique_fit_method_descriptors[
  !(unique_fit_method_descriptors %in% unique_coef_name_component_descriptors)]
intercept_obs_dist_complex <- unique_obs_dist_complex_descriptors[
  !(unique_obs_dist_complex_descriptors %in% unique_coef_name_component_descriptors)]
intercept_time_dep <- unique_time_dep_descriptors[
  !(unique_time_dep_descriptors %in% unique_coef_name_component_descriptors)]
for(coef_ind in seq(from = 1, to = length(coef_names))) {
	split_name <- unlist(strsplit(coef_names[[coef_ind]], ":"))
	if(!any(split_name %in% unique_fit_method_descriptors[unique_fit_method_descriptors != intercept_fit_method])) {
		split_name <- c(split_name, unique_fit_method_descriptors)
	}
	if(!any(split_name %in% unique_obs_dist_complex_descriptors[unique_obs_dist_complex_descriptors != intercept_obs_dist_complex])) {
		split_name <- c(split_name, unique_obs_dist_complex_descriptors)
	}
	if(!any(split_name %in% unique_time_dep_descriptors[unique_time_dep_descriptors != intercept_time_dep])) {
		split_name <- c(split_name, unique_time_dep_descriptors)
	}

	lc_df[[paste0("coef", coef_ind)]] <- 0
	lc_df[[paste0("coef", coef_ind)]][
	  lc_df$fit_method_descriptor %in% split_name &
		lc_df$obs_dist_complex_descriptor %in% split_name &
	  lc_df$time_dep_descriptor %in% split_name] <- 1
}

## contrasts of
## (mean performance method 1) - (mean performance method 2) for all pairs of methods
## within obs_dist_complex and time_dep
rowind <- nrow(lc_df) # index of new row to add to lc_df
confint_rows <- c() # rows for which to compute confidence intervals
for(fit_method1_ind in seq(from = 1, to = length(unique_fit_methods) - 1)) {
  for(fit_method2_ind in seq(from = fit_method1_ind + 1, to = length(unique_fit_methods))) {
    fit_method1 <- unique_fit_methods[fit_method1_ind]
    fit_method2 <- unique_fit_methods[fit_method2_ind]
    
    for(obs_dist_complex_val in unique_obs_dist_complex) {
      for(time_dep_val in unique_time_dep) {
        rowind <- rowind + 1
      	confint_rows <- c(confint_rows, rowind)
    	  
        m1_rowind <- which(lc_df$name == paste0(fit_method1, "-", obs_dist_complex_val, "-", time_dep_val, "-fit_method", fit_method1))
        m2_rowind <- which(lc_df$name == paste0(fit_method2, "-", obs_dist_complex_val, "-", time_dep_val, "-fit_method", fit_method2))
        
      	lc_df[rowind, ] <- rep(NA, ncol(lc_df))
        lc_df$name[rowind] <- paste0(fit_method1, "-", fit_method2, "-", obs_dist_complex_val, "-", time_dep_val)
      	lc_df$obs_dist_complex[rowind] <- obs_dist_complex_val
      	lc_df$time_dep[rowind] <- time_dep_val
      	lc_df[rowind, coef_cols] <- lc_df[m1_rowind, coef_cols] - lc_df[m2_rowind, coef_cols]
      }
    }
  }
}

lc_df$name <- factor(lc_df$name, levels = lc_df$name)

K_mat <- as.matrix(lc_df[, coef_cols])

# get point estimates
lc_df$pt_est <- as.vector(K_mat %*% matrix(fixef(results_fit)))

# get familywise CIs
lc_df$fam_CI_lb <- NA
lc_df$fam_CI_ub <- NA
fam_CI_obj <- glht(results_fit, linfct = K_mat[confint_rows, ])
temp <- confint(fam_CI_obj)$confint
lc_df$fam_CI_lb[confint_rows] <- temp[, 2]
lc_df$fam_CI_ub[confint_rows] <- temp[, 3]

# get individual CIs
lc_df$ind_CI_lb <- NA
lc_df$ind_CI_ub <- NA
for(rowind in confint_rows) {
	ind_CI_obj <- glht(results_fit, linfct = K_mat[rowind, , drop = FALSE])
	temp <- confint(ind_CI_obj)$confint
	lc_df$ind_CI_lb[rowind] <- temp[, 2]
	lc_df$ind_CI_ub[rowind] <- temp[, 3]
}


summary_figure_df <-
  lc_df[21:60, c("name", "obs_dist_complex", "time_dep", "pt_est", "fam_CI_lb", "fam_CI_ub", "ind_CI_lb", "ind_CI_ub")]
summary_figure_df$method_contrast <- 
  sapply(strsplit(as.character(summary_figure_df$name), "-", fixed = TRUE), function(comp) { paste(comp[1], "-", comp[2]) })

ggplot(data = summary_figure_df) +
  geom_point(aes(x = method_contrast, y = pt_est)) +
  geom_errorbar(aes(x = method_contrast, ymin = fam_CI_lb, ymax = fam_CI_ub)) +
  facet_grid(time_dep ~ obs_dist_complex) +
  xlab("Classification Model Pair") +
  ylab("Difference in Proportion Correct") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
