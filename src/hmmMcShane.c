/* 
 ============================================================================
 Name        : hmmMcShane.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : Implementation of McShane's discriminative HMM formulation
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

#include <math.h>

#include "utility.h"
#include "crftype.h"
#include "lccrf.h"
#include "hmmMcShane.h"



void calc_log_McShane_class_probs_given_log_static_class_probs(CRF *HMM_params, int N, int *T, double *log_marginal_class_probs,
		double **forward_vars, double **backward_vars, double **log_prop_obs_prob_given_state, double ***log_static_class_probs, double ***log_mcshane_class_probs) {
	for(int i = 0; i < N; i++) {
		Rprintf("prop log obs probs:\n");
		for(int t = 0; t < *(T + i); t++) {
			Rprintf("t = %i: ", t);
			for(int s = 0; s < HMM_params->S; s++) {
				*(*(log_prop_obs_prob_given_state + t) + s) = *(*(*(log_static_class_probs + i) + t) + s) - *(log_marginal_class_probs + s);
				Rprintf("%lf, ", *(*(log_prop_obs_prob_given_state + t) + s));
			}
			Rprintf("\n");
		}

		calc_forward_vars_logs(HMM_params, *(T + i), log_prop_obs_prob_given_state, forward_vars);
		calc_backward_vars_logs(HMM_params, *(T + i), log_prop_obs_prob_given_state, backward_vars);

		Rprintf("\nforward vars:\n");
		for(int t = 0; t < *(T + i); t++) {
			Rprintf("t = %i: ", t);
			for(int s = 0; s < HMM_params->S; s++) {
				Rprintf("%lf, ", *(*(forward_vars + t) + s));
			}
			Rprintf("\n");
		}

		Rprintf("\nbackward vars:\n");
		for(int t = 0; t < *(T + i); t++) {
			Rprintf("t = %i: ", t);
			for(int s = 0; s < HMM_params->S; s++) {
				Rprintf("%lf, ", *(*(backward_vars + t) + s));
			}
			Rprintf("\n");
		}


		// calculate the probabilities of each class at all t based on all models combined
		calc_ps_logs(HMM_params, *(T + i), forward_vars, backward_vars, *(log_mcshane_class_probs + i));
	}
}



SEXP calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(
		SEXP HMM_Xptr,
		SEXP log_marginal_class_probsp,
		SEXP Np,
		SEXP Tp,
		SEXP log_static_class_probsp) {
	// DEFINE SOME INDEX AND TEMPORARY VARIABLES
	int i, t, s, maxT, nRprotect = 0, N = *INTEGER(Np);
	int *T = INTEGER(Tp);
	double *tempdblptr;
	

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// HMM_params contains the parameters for the HMM.
	CRF *HMM_params = CRF_from_RXptr(HMM_Xptr);

	// a variable to contain the return value
	SEXP retval;
	retval = PROTECT(allocVector(VECSXP, N));
	for(i = 0; i < N; i++) {
		SET_VECTOR_ELT(retval, i, allocMatrix(REALSXP, *(T + i), HMM_params->S));
	}
	nRprotect++;


	
	// prev_fs stores the combined fs from previous additive components.
	// note that we want to copy the values to new memory instead of pointing to the R objects
	// because we pass prev_fs to IACL_log_lik_one_sub by reference so that it can modify it
	// in case we need to return the combined fs from all models
	double ***log_static_class_probs = Calloc(N, double**);
	double ***log_mcshane_class_probs = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(log_static_class_probs + i) = Calloc( *(T + i), double* );
		*(log_mcshane_class_probs + i) = Calloc( *(T + i), double* );
		for(t = 0; t < T[i]; t++) {
			*(*(log_static_class_probs + i) + t) = Calloc(HMM_params->S, double);
			*(*(log_mcshane_class_probs + i) + t) = Calloc(HMM_params->S, double);
		}
		
		tempdblptr = REAL(VECTOR_ELT(log_static_class_probsp, i));
		for(t = 0; t < *(T + i); t++) {
			for(s = 0; s < HMM_params->S; s++) {
				*(*(*(log_static_class_probs + i) + t) + s) = *(tempdblptr + t + s * *(T + i));
			}
		}
	}

	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = *(T);
	for(i = 1; i < N; i++)
		if(maxT < *(T + i))
			maxT = *(T + i);

	double **forward_vars = Calloc(maxT, double*);
	double **backward_vars = Calloc(maxT, double*);
	double **obs_prob_given_state = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(HMM_params->S, double);
		backward_vars[t] = Calloc(HMM_params->S, double);
		obs_prob_given_state[t] = Calloc(HMM_params->S, double);
	}
	
	
	// DO THE CALCULATION OF THE CLASS PROBABILITIES
	calc_log_McShane_class_probs_given_log_static_class_probs(HMM_params, N, T, REAL(log_marginal_class_probsp),
			forward_vars, backward_vars, obs_prob_given_state, log_static_class_probs, log_mcshane_class_probs);
	
	
	for(i = 0; i < N; i++) {
		tempdblptr = REAL(VECTOR_ELT(retval, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < HMM_params->S; s++) {
				*(tempdblptr + t + s * *(T + i)) = *(*(*(log_mcshane_class_probs + i) + t) + s);
			}
		}
	}

	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
		Free(backward_vars[t]);
		Free(obs_prob_given_state[t]);
	}
	Free(forward_vars);
	Free(backward_vars);
	Free(obs_prob_given_state);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(log_mcshane_class_probs + i) + t));
			Free(*(*(log_static_class_probs + i) + t));
		}
		Free(*(log_mcshane_class_probs + i));
		Free(*(log_static_class_probs + i));
	}
	Free(log_mcshane_class_probs);
	Free(log_static_class_probs);

	UNPROTECT(nRprotect);

	return retval;
}
