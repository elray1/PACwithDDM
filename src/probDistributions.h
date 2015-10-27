/*
 ============================================================================
 Name        : probDistributions.h
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : probability distributions
 ============================================================================
 */

double ddirich(double *x, double *alpha, int len);

void rdirich(double *alpha, int len, double *target);

double dMVN(double *x, int length, int x_pointer_offset, double *mu, double *Sigma, double *log_det_Sigma, int retlog);

SEXP dMVN_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog);

void dMVN_multiobs(double *x, int n, int d, double *mu, double *Sigma, double *log_det_Sigma, double *result, int retlog);

SEXP dMVN_Prec_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog);

void dMVN_Prec_multiobs(double *x, int n, int d, double *mu, double *Prec, double *log_det_Sigma, double *result, int retlog);

SEXP dMVN_DiagPrec_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog);

void dMVN_DiagPrec_multiobs(double *x, int n, int d, double *mu, double *Prec, double *log_det_Sigma, double *result, int retlog);

void rMVN(int length, double *mu, double *Sigma, double *target);

SEXP dGMM_Prec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precsp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog);

SEXP dGMM_samePrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precp, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog);

SEXP dGMM_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Sigmasp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog);

SEXP dGMM_sameSigma_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Sigmap, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog);

SEXP dGMM_DiagPrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precsp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog);

SEXP dGMM_sameDiagPrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precp, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog);

int uneq_prob_sample(double *prop_prob_vec, int length);

int runif_int(int min, int max);

void sample_int_without_replace_FisherYates(int *values, int K, int *result, int N);

void sample_int_without_replace(int min, int max, int N, int *result);
