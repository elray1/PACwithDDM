#ifndef LCCRF_H
#define LCCRF_H

#include "crftype.h"

void calc_forward_vars_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **forward_vars);

void calc_backward_vars_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **backward_vars);

void calc_ps_logs(CRF *CRF_params, int T, double **forward_vars, double **backward_vars, double **ps);

void calc_jtps_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **forward_vars, double **backward_vars, double **jtps);

void CRF_log_lik_given_obs_probs(CRF *CRF_params, int N, int *T, double ***log_obs_prob_given_state, int **a,
		double *loglik, double **forward_vars);

SEXP marginal_CRF_log_lik_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP log_obs_probs_by_statep, SEXP retallp);

SEXP CRF_log_lik_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP log_obs_probs_by_statep);

double partial_CRF_log_lik_wrt_omega_r_s(CRF *CRF_params, int N, int *T, int **a, double ***jtps, int r, int s);

double partial_CRF_log_lik_wrt_omega_reduced_parameterization(CRF *CRF_params, int N, int *T, int **a, double ***jtps);

SEXP CRF_gradient_log_lik_wrt_omegas_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP combined_log_obs_probs_by_statep, SEXP reduced_trans_mat_parameterizationp);

void calc_log_obs_prob_one_sub_given_state_globalSeqCRF(CRF *CRF_params, int T, double *betas, int *cols, int C,
		double *X, double **log_obs_prob_given_state);

SEXP globalSeqCRF_update_log_lik_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP prev_log_obs_probs_by_statep,
		SEXP betasp, SEXP colsp, SEXP retallp);

SEXP globalSeqCRF_update_gradient_log_lik_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP prev_log_obs_probs_by_statep,
		SEXP betasp, SEXP colsp, SEXP beta_diff_colsp, SEXP diff_wrt_omegasp, SEXP reduced_log_mat_parameterizationp);

#endif
