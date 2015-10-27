#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "lccrf.h"
#include "hmmMcShane.h"
#include "crftype.h"
#include "utility.h"
#include "signedlog.h"
#include "probDistributions.h"


////////////////////////////////////////////////////////////////////////////////
// Basic calculations for linear chain crfs/hmms
////////////////////////////////////////////////////////////////////////////////

// Calculate logarithm of forward variables at all time points for a given subject
//
// INPUTS:
// CRF_params is an object with all CRF parameters
// T is the number of time points for the given subject
// prob_obs_t_given_state_s is a pointer to an array of length T_i, where entry i is a S by 1 matrix
//      containing the probability of the data given state 1, ..., S
//
// OUTPUTS:
// No return value.
// The forward_vars argument is modified to contain the logs of the calculated forward vars.
void calc_forward_vars_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **forward_vars) {
	int t, r, s;

	// initialize
	for(s = 0; s < CRF_params->S; s++)
		*(*(forward_vars + 0) + s) = *(CRF_params->log_pi + s) + *(*(prob_obs_t_given_state_s + 0) + s);

	for(t = 1; t < T; t++) {
		for(s = 0; s < CRF_params->S; s++) {
			*(*(forward_vars + t) + s) = *(*(forward_vars + t - 1) + 0) + *(CRF_params->log_trans_matrix + 0 + s * CRF_params->S) + *(*(prob_obs_t_given_state_s + t) + s);
			for(r = 1; r < CRF_params->S; r++) {
				*(*(forward_vars + t) + s) = logspace_add_safe(*(*(forward_vars + t) + s), *(*(forward_vars + t - 1) + r) + *(CRF_params->log_trans_matrix + r + s * CRF_params->S) + *(*(prob_obs_t_given_state_s + t) + s));
			}
		}
	}
}



// Calculate logarithm of backward variables at all time points for a given subject
//
// INPUTS:
// CRF_params is an object with all CRF parameters
// T is the number of time points for the given subject
// prob_obs_t_given_state_s is a pointer to an array of length T_i, where entry i is a S by 1 matrix
//      containing the probability of the data given state 1, ..., S
//
// OUTPUTS:
// No return value.
// The backward_vars argument is modified to contain the logs of the calculated backward vars.
void calc_backward_vars_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **backward_vars) {
	int t, r, s;

	// initialize
	t = T - 1; // T_sub - 1
	for(s = 0; s < CRF_params->S; s++) {
		*(*(backward_vars + t) + s) = 0;
	}

	// recursion
	for(t = T - 2; t >= 0 ; t--) {
		for(r = 0; r < CRF_params->S; r++) {
			*(*(backward_vars + t) + r) = *(*(backward_vars + t + 1) + 0) + *(CRF_params->log_trans_matrix + r + 0 * CRF_params->S) + *(*(prob_obs_t_given_state_s + t + 1) + 0);

			for(s = 1; s < CRF_params->S; s++) {
				*(*(backward_vars + t) + r) = logspace_add_safe(*(*(backward_vars + t) + r), *(*(backward_vars + t + 1) + s) + *(CRF_params->log_trans_matrix + r + s * CRF_params->S) + *(*(prob_obs_t_given_state_s + t + 1) + s));
			}
		}
	}
}



// Calculate logarithm of the class probabilities at each time point (marginalized over all other time points)
// 
// INPUTS:
// CRF_params is an object with all CRF parameters
// T is the number of time points for the given subject
// forward_vars, backward_vars are arrays containing the logarithms of the forward and backward variables
// 
// OUTPUTS:
// No return value.
// ps is modified to contain the marginal class probabilies at each time point
void calc_ps_logs(CRF *CRF_params, int T, double **forward_vars, double **backward_vars, double **ps) {
	int t, s;

	// normalization constant
	double tempdbl;
	tempdbl = *(*(forward_vars + T - 1) + 0);
	for(s = 1; s < CRF_params->S; s++)
		tempdbl = logspace_add_safe(tempdbl, *(*(forward_vars + T - 1) + s));

	// marginal class probabilities at each time point
	for(t = 0; t < T; t++) {
		for(s = 0; s < CRF_params->S; s++)
			*(*(ps + t) + s) = *(*(forward_vars + t) + s) + *(*(backward_vars + t) + s) - tempdbl;

	}
}



// Calculate logarithms of joint probabilities of class r at time t-1 and s at time t for each time point t and
// combination of classes r and s, marginalizing over classes for all other time points.
//
// INPUTS:
// CRF_params is an object with all CRF parameters
// T is the number of time points for the given subject
// forward_vars, backward_vars are arrays containing the logarithms of the forward and backward variables
//
// OUTPUTS:
// No return value.
// jtps is modified to contain the logarithm of the marginal class probabilies of
// state r at time t-1 and s at time t
void calc_jtps_logs(CRF *CRF_params, int T, double **prob_obs_t_given_state_s, double **forward_vars, double **backward_vars, double **jtps) {
	int t, r, s;

	// normalization constant
	double tempdbl;
	tempdbl = *(*(forward_vars + T - 1) + 0);
	for(s = 1; s < CRF_params->S; s++)
		tempdbl = logspace_add_safe(tempdbl, *(*(forward_vars + T - 1) + s));

	// for each time t = 1, ..., T-1 (recalling we are indexing from 0),
	// log of probability state r at t - 1, s at t
	for(t = 1; t < T; t++) {
		for(r = 0; r < CRF_params->S; r++)
			for(s = 0; s < CRF_params->S; s++)
				*(*(jtps + t) + r * CRF_params->S + s) = *(*(forward_vars + t - 1) + r) + *(CRF_params->log_trans_matrix + r + s * CRF_params->S) + *(*(prob_obs_t_given_state_s + t) + s) + *(*(backward_vars + t) + s) - tempdbl;
	}
}



// Calculate logarithm of conditional probability of a sequence of states for one subject given observations.
//
// @param CRF_params is an object with all CRF parameters
// @param T is the number of time points for the given subject
// @param log_obs_probs_given_state 
// @param a a vector of observed states with length T
// @param forward_vars an array containing the logarithms of the forward variables
//
// @return log(conditional probability of state sequence a given data)
double CRF_sequence_log_prob(CRF *CRF_params, int T, double **log_obs_probs_given_state, int *a, double **forward_vars) {
	double result, Z;

	// initialize
	result = *(CRF_params->log_pi + *(a + 0) - 1) + *(*(log_obs_probs_given_state + 0) + *(a + 0) - 1);

	// update
	for(int t = 1; t < T; t++)
		result += *(CRF_params->log_trans_matrix + *(a + t - 1) - 1 + (*(a + t) - 1) * CRF_params->S) + *(*(log_obs_probs_given_state + t) + *(a + t) - 1);

	// normalizing factor
	Z = *(*(forward_vars + T - 1));
	for(int s = 1; s < CRF_params->S; s++)
		Z = logspace_add_safe(*(*(forward_vars + T - 1) + s), Z);

	result = result - Z;

	return result;
}



/**
 * Log-likelihood calculation for one subject's data
 * 
 * @param CRF_params
 * @param T
 * @param log_obs_probs_given_state
 * @param a
 * @param log_lik
 * @param forward_vars
 * 
 * @return none 
 */
void CRF_log_lik_given_obs_probs_one_subject(CRF *CRF_params, int T, double **log_obs_probs_given_state, int *a, double *loglik, double **forward_vars) {
	// calculate the forward/backward vars
	calc_forward_vars_logs(CRF_params, T, log_obs_probs_given_state, forward_vars);

	// update loglikelihood with contribution from this subject
	*loglik += CRF_sequence_log_prob(CRF_params, T, log_obs_probs_given_state, a, forward_vars);
}



/**
 * Log-likelihood calculation for multiple subjects' data
 * 
 * @param CRF_params
 * @param N
 * @param T
 * @param log_obs_probs_given_state
 * @param a
 * @param log_lik
 * @param forward_vars
 * 
 * @return none
 */
void CRF_log_lik_given_obs_probs(CRF *CRF_params, int N, int *T, double ***log_obs_prob_given_state, int **a,
		double *loglik, double **forward_vars) {
	int i;
	
	*loglik = 0;

	for(i = 0; i < N; i++)
		CRF_log_lik_given_obs_probs_one_subject(CRF_params, *(T + i), *(log_obs_prob_given_state + i), *(a + i),
			loglik, forward_vars);
}





////////////////////////////////////////////////////////////////////////////////
// Functions to do set-up specific to parametric CRFs
////////////////////////////////////////////////////////////////////////////////


// Calculate the log probabilities of the observed data given the states at all times and each state for one subject
// based on a parametric linear chain CRF specification with given parameter values
// 
// INPUTS:
// CRF_params is an object with all CRF parameters
// T is the number of time points for the given subject
// betas vector of non-zero beta values associated with each state, arranged as a D by C matrix
// cols columns of observed data X with non-zero beta coefficients
// C length of cols
// X is a pointer to a T by m matrix, where row t contains the observed values of the m covariates
//      x_{t, 1}, ..., x_{t, m} for a particular subject at time t
// log_obs_prob_given_state is a pointer to an array of length T_i, where entry i is a S by 1 matrix
//			that is not yet populated
// 
// OUTPUTS:
// No return value.
// log_obs_prob_given_state is modified to contain
//      the probability of the data X_t given state s for t = 1, ..., T_i and s = 1, ..., S
void calc_log_obs_prob_one_sub_given_state_parametric_CRF(CRF *CRF_params, int T, double *betas, int *cols, int C,
		double *X, double **log_obs_prob_given_state) {
	int t, s, c;
	double raw_log_obs_prob;

	// at each state, compute the likelihood for the observed features at time t given the parameters for that state
	for(t = 0; t < T; t++) {
		for(s = 0; s < CRF_params->S; s++) {
			raw_log_obs_prob = 0;
			for(c = 0; c < C; c++) {
				raw_log_obs_prob += *(betas + s * CRF_params->D + *(cols + c) - 1) * *(X + t + (*(cols + c) - 1) * T);
			}
			
			*(*(log_obs_prob_given_state + t) + s) = raw_log_obs_prob;
		}
	}
}



void update_log_obs_probs_all_subs_given_state_parametric_CRF(CRF *CRF_params, int N, int *T, double *betas, int *cols, int C,
		double **X, double ***prev_log_obs_prob_given_state, double ***log_obs_prob_given_state) {
	int i, t, s;
	for(i = 0; i < N; i++) {
		calc_log_obs_prob_one_sub_given_state_parametric_CRF(CRF_params, *(T + i), betas, cols, C, *(X + i), *(log_obs_prob_given_state + i));

		for(t = 0; t < *(T + i); t++)
			for(s = 0; s < CRF_params->S; s++)
				*(*(*(log_obs_prob_given_state + i) + t) + s) = *(*(*(log_obs_prob_given_state + i) + t) + s) + *(*(*(prev_log_obs_prob_given_state + i) + t) + s);
	}
}



////////////////////////////////////////////////////////////////////////////////
// R interfaces for log-likelihood calculation functions
////////////////////////////////////////////////////////////////////////////////

SEXP parametric_CRF_update_log_lik_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP prev_log_obs_probs_by_statep,
		SEXP betasp, SEXP colsp, SEXP retallp) {
	// DEFINE SOME INDEX AND TEMPORARY VARIABLES
	int i, t, s, maxT, nRprotect = 0, N = length(observedCRFsp), C = length(colsp), retall = *INTEGER(retallp);
	int T[N];
	double *tempdblptr1;

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// CRF_params contains the parameters for the CRF.
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);

	for(i = 0; i < N; i++) {
		T[i] = *(INTEGER(getAttrib(getListElement(VECTOR_ELT(observedCRFsp, i), "X"), R_DimSymbol)));
	}

	// a variable to contain the return value
	SEXP retval;
	if(retall) {
		retval = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(retval, 0, allocVector(REALSXP, 1));
		SET_VECTOR_ELT(retval, 1, allocVector(VECSXP, N));
		for(i = 0; i < N; i++) {
			SET_VECTOR_ELT(VECTOR_ELT(retval, 1), i, allocMatrix(REALSXP, T[i], CRF_params->S));
		}
		nRprotect++;
	} else {
		retval = PROTECT(allocVector(REALSXP, 1));
		nRprotect++;
	}


	
	double ***prev_log_obs_probs_by_state = Calloc(N, double**);
	double ***log_obs_probs_by_state = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(prev_log_obs_probs_by_state + i) = Calloc( T[i], double* );
		*(log_obs_probs_by_state + i) = Calloc( T[i], double* );
		for(t = 0; t < T[i]; t++) {
			*(*(prev_log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
			*(*(log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
		}
		
		tempdblptr1 = REAL(VECTOR_ELT(prev_log_obs_probs_by_statep, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < CRF_params->S; s++) {
				*(*(*(prev_log_obs_probs_by_state + i) + t) + s) = *(tempdblptr1 + t + s * T[i]);
			}
		}
	}

	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = T[0];
	for(i = 1; i < N; i++)
		if(maxT < T[i])
			maxT = T[i];

	double **forward_vars = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(CRF_params->S, double);
	}
	
	// pointers to the observed class and covariates at all time points for each subject
	// the observed class
	int *a[N];
	double *X[N];

	for(i = 0; i < N; i++) {
		a[i] = INTEGER(getListElement(VECTOR_ELT(observedCRFsp, i), "y"));
		X[i] = REAL(getListElement(VECTOR_ELT(observedCRFsp, i), "X"));
	}
	
	// DO THE CALCULATION OF THE LOG-LIKELIHOOD
	// calculate updated log_obs_probs
	update_log_obs_probs_all_subs_given_state_parametric_CRF(CRF_params, N, T, REAL(betasp), INTEGER(colsp), C,
		X, prev_log_obs_probs_by_state, log_obs_probs_by_state);

	// the calculated log-likelihood
	double loglik;
	
	CRF_log_lik_given_obs_probs(CRF_params, N, T, log_obs_probs_by_state, a,
		&loglik, forward_vars);
	
	if(retall) {
		// if requested, store the combined fs in the return value.
		for(i = 0; i < N; i++) {
			tempdblptr1 = REAL(VECTOR_ELT(VECTOR_ELT(retval, 1), i));
			for(t = 0; t < T[i]; t++) {
				for(s = 0; s < CRF_params->S; s++) {
					*(tempdblptr1 + t + s * T[i]) = *(*(*(log_obs_probs_by_state + i) + t) + s);
				}
			}
		}

		*REAL(VECTOR_ELT(retval,0)) = loglik;
	} else {
		*REAL(retval) = loglik;
	}


	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
	}
	Free(forward_vars);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(prev_log_obs_probs_by_state + i) + t));
			Free(*(*(log_obs_probs_by_state + i) + t));
		}
		Free(*(prev_log_obs_probs_by_state + i));
		Free(*(log_obs_probs_by_state + i));
	}
	Free(prev_log_obs_probs_by_state);
	Free(log_obs_probs_by_state);

	UNPROTECT(nRprotect);

	return retval;
}





SEXP CRF_log_lik_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP log_obs_probs_by_statep) {
	// DEFINE SOME INDEX AND TEMPORARY VARIABLES
	int i, t, s, maxT, nRprotect = 0, N = length(observedCRFsp);
	int T[N];
	double *tempdblptr1;

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// CRF_params contains the parameters for the CRF.
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);

	for(i = 0; i < N; i++) {
		T[i] = *(INTEGER(getAttrib(getListElement(VECTOR_ELT(observedCRFsp, i), "X"), R_DimSymbol)));
	}

	// a variable to contain the return value
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	nRprotect++;


	
	double ***log_obs_probs_by_state = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(log_obs_probs_by_state + i) = Calloc( T[i], double* );
		for(t = 0; t < T[i]; t++) {
			*(*(log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
		}
		
		tempdblptr1 = REAL(VECTOR_ELT(log_obs_probs_by_statep, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < CRF_params->S; s++) {
				*(*(*(log_obs_probs_by_state + i) + t) + s) = *(tempdblptr1 + t + s * T[i]);
			}
		}
	}

	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = T[0];
	for(i = 1; i < N; i++)
		if(maxT < T[i])
			maxT = T[i];

	double **forward_vars = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(CRF_params->S, double);
	}
	
	// pointers to the observed class and covariates at all time points for each subject
	// the observed class
	int *a[N];

	for(i = 0; i < N; i++)
		a[i] = INTEGER(getListElement(VECTOR_ELT(observedCRFsp, i), "y"));
	
	// DO THE CALCULATION OF THE LOG-LIKELIHOOD
	// the calculated log-likelihood
	double *loglik;
	loglik = REAL(retval);
	
	CRF_log_lik_given_obs_probs(CRF_params, N, T, log_obs_probs_by_state, a,
		loglik, forward_vars);
	

	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
	}
	Free(forward_vars);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(log_obs_probs_by_state + i) + t));
		}
		Free(*(log_obs_probs_by_state + i));
	}
	Free(log_obs_probs_by_state);

	UNPROTECT(nRprotect);

	return retval;
}



void marginal_CRF_log_lik_given_obs_probs_one_subject(CRF *CRF_params, int T, double **log_obs_probs_given_state, int *a, double *loglik,
		double **forward_vars, double **backward_vars, double **ps) {
	// VARIABLE DECALARATIONS

	// temporary and index variables
	int t;

	// DO THE CALCULATION OF THE LIKELIHOOD

	// calculate the forward/backward vars
	calc_forward_vars_logs(CRF_params, T, log_obs_probs_given_state, forward_vars);
	calc_backward_vars_logs(CRF_params, T, log_obs_probs_given_state, backward_vars);

	// calculate the probabilities of each class at all t based on all models combined
	calc_ps_logs(CRF_params, T, forward_vars, backward_vars, ps);

	for(t = 0; t < T; t++)
		*loglik += *(*(ps + t) + *(a + t) - 1);
}

void marginal_CRF_log_lik_given_obs_probs(CRF *CRF_params, int N, int *T, double ***log_obs_prob_given_state, int **a,
		double *loglik, double **forward_vars, double **backward_vars, double ***ps) {
	int i;
	
	*loglik = 0;

	for(i = 0; i < N; i++)
		marginal_CRF_log_lik_given_obs_probs_one_subject(CRF_params, *(T + i), *(log_obs_prob_given_state + i), *(a + i),
			loglik, forward_vars, backward_vars, *(ps + i));
}

SEXP marginal_CRF_log_lik_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP log_obs_probs_by_statep, SEXP retallp) {
	// DEFINE SOME INDEX AND TEMPORARY VARIABLES
	int i, t, s, maxT, nRprotect = 0, N = length(observedCRFsp);
	int T[N];
	double *tempdblptr1;

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// CRF_params contains the parameters for the CRF.
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);

	for(i = 0; i < N; i++) {
		T[i] = *(INTEGER(getAttrib(getListElement(VECTOR_ELT(observedCRFsp, i), "X"), R_DimSymbol)));
	}

	// a variable to contain the return value
	int retall = *INTEGER(retallp);
	SEXP retval;
	if(retall) {
		retval = PROTECT(allocVector(VECSXP, 2));
		nRprotect++;
		SET_VECTOR_ELT(retval, 0, allocVector(REALSXP, 1));
		SET_VECTOR_ELT(retval, 1, allocVector(VECSXP, N));
		for(i = 0; i < N; i++) {
			SET_VECTOR_ELT(VECTOR_ELT(retval, 1), i, allocMatrix(REALSXP, T[i], CRF_params->S));
			tempdblptr1 = REAL(VECTOR_ELT(VECTOR_ELT(retval, 1), i));
		}
	} else {
		retval = PROTECT(allocVector(REALSXP, 1));
		nRprotect++;
	}

	
	double ***log_obs_probs_by_state = Calloc(N, double**);
	double ***ps = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(log_obs_probs_by_state + i) = Calloc( T[i], double* );
		*(ps + i) = Calloc( T[i], double* );
		for(t = 0; t < T[i]; t++) {
			*(*(log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
			*(*(ps + i) + t) = Calloc(CRF_params->S, double);
		}
		
		tempdblptr1 = REAL(VECTOR_ELT(log_obs_probs_by_statep, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < CRF_params->S; s++) {
				*(*(*(log_obs_probs_by_state + i) + t) + s) = *(tempdblptr1 + t + s * T[i]);
			}
		}
	}

	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = T[0];
	for(i = 1; i < N; i++)
		if(maxT < T[i])
			maxT = T[i];

	double **forward_vars = Calloc(maxT, double*);
	double **backward_vars = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(CRF_params->S, double);
		backward_vars[t] = Calloc(CRF_params->S, double);
	}
	
	// pointers to the observed class and covariates at all time points for each subject
	// the observed class
	int *a[N];

	for(i = 0; i < N; i++)
		a[i] = INTEGER(getListElement(VECTOR_ELT(observedCRFsp, i), "y"));
	
	// DO THE CALCULATION OF THE LOG-LIKELIHOOD
	// the calculated log-likelihood
	double *loglik;
	if(retall) {
		loglik = REAL(VECTOR_ELT(retval, 0));
	} else {
		loglik = REAL(retval);
	}
	
	marginal_CRF_log_lik_given_obs_probs(CRF_params, N, T, log_obs_probs_by_state, a,
		loglik, forward_vars, backward_vars, ps);
	
	if(retall) {
		// if requested, store the log class probabilities in the return value.
		for(i = 0; i < N; i++) {
			tempdblptr1 = REAL(VECTOR_ELT(VECTOR_ELT(retval, 1), i));
			for(t = 0; t < T[i]; t++) {
				for(s = 0; s < CRF_params->S; s++) {
					*(tempdblptr1 + t + s * T[i]) = *(*(*(ps + i) + t) + s);
				}
			}
		}
	}

	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
		Free(backward_vars[t]);
	}
	Free(forward_vars);
	Free(backward_vars);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(log_obs_probs_by_state + i) + t));
			Free(*(*(ps + i) + t));
		}
		Free(*(log_obs_probs_by_state + i));
		Free(*(ps + i));
	}
	Free(log_obs_probs_by_state);
	Free(ps);

	UNPROTECT(nRprotect);

	return retval;
}





////////////////////////////////////////////////////////////////////////////////
// Gradient calculations
////////////////////////////////////////////////////////////////////////////////

double partial_parametric_CRF_log_lik_wrt_beta_s_d_given_obs_probs(CRF *CRF_params, int N, int *T, int **a, double **X, double ***ps, int s, int d) {
	signedlog *resultsl = new_signedlog(-INFINITY, 1), *newsl = new_signedlog(-INFINITY, 1);
	double temp;
	int T_i;

	for(int i = 0; i < N; i++) {
		T_i = *(T + i);
		for(int t = 0; t < T_i; t++) {
			if(*(*(a + i) + t) - 1 == s) {
				newsl->slval = logspace_sub(0.0, *(*(*(ps + i) + t) + s));
				if(isnan(newsl->slval))
					newsl->slval = -INFINITY;
				newsl->slsign = 1;
			} else {
				newsl->slval = *(*(*(ps + i) + t) + s);
				newsl->slsign = -1;
			}

			logspace_mult_signedlogs(newsl->slval, newsl->slsign, log(fabs(*(*(X + i) + t + d * T_i))), (int)copysignf(1.0, *(*(X + i) + t + d * T_i)), &(newsl->slval), &(newsl->slsign));

			logspace_add_signedlogs(resultsl->slval, resultsl->slsign, newsl->slval, newsl->slsign, &(resultsl->slval), &(resultsl->slsign));
		}
	}

	temp = signedlog_to_dbl(resultsl);

	free_signedlog(resultsl);
	free_signedlog(newsl);

	return(temp);
}

SEXP parametric_CRF_update_gradient_log_lik_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP prev_log_obs_probs_by_statep,
		SEXP betasp, SEXP colsp, SEXP beta_diff_colsp, SEXP diff_wrt_omegasp, SEXP reduced_log_mat_parameterizationp) {
	// DEFINE SOME INDEX AND TEMPORARY VARIABLES
	int i, t, r, s, d, maxT, nRprotect = 0, N = length(observedCRFsp), C = length(colsp), num_betas_diff = length(beta_diff_colsp),
		diff_wrt_omegas = *INTEGER(diff_wrt_omegasp), reduced_log_mat_parameterization = *INTEGER(reduced_log_mat_parameterizationp);
	int *beta_diff_cols = INTEGER(beta_diff_colsp);
	int T[N];
	double *tempdblptr1;

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// CRF_params contains the parameters for the CRF.
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);
	int Sm1 = CRF_params->S - 1;

	for(i = 0; i < N; i++) {
		T[i] = *(INTEGER(getAttrib(getListElement(VECTOR_ELT(observedCRFsp, i), "X"), R_DimSymbol)));
	}

	// a variable to contain the return value
	SEXP retval;
	if(diff_wrt_omegas) {
		if(reduced_log_mat_parameterization) {
			retval = PROTECT(allocVector(REALSXP, 1 + Sm1 * num_betas_diff));
		} else {
			retval = PROTECT(allocVector(REALSXP, CRF_params->S * CRF_params->S - 1 + Sm1 * num_betas_diff));
		}
	} else {
		retval = PROTECT(allocVector(REALSXP, Sm1 * num_betas_diff));
	}
	nRprotect++;


	
	double ***prev_log_obs_probs_by_state = Calloc(N, double**);
	double ***log_obs_probs_by_state = Calloc(N, double**);
	double ***ps = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(prev_log_obs_probs_by_state + i) = Calloc( T[i], double* );
		*(log_obs_probs_by_state + i) = Calloc( T[i], double* );
		*(ps + i) = Calloc( T[i], double* );
		for(t = 0; t < T[i]; t++) {
			*(*(prev_log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
			*(*(log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
			*(*(ps + i) + t) = Calloc(CRF_params->S, double);
		}
		
		tempdblptr1 = REAL(VECTOR_ELT(prev_log_obs_probs_by_statep, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < CRF_params->S; s++) {
				*(*(*(prev_log_obs_probs_by_state + i) + t) + s) = *(tempdblptr1 + t + s * T[i]);
			}
		}
	}

	double ***jtps;
	if(diff_wrt_omegas) {
		jtps = Calloc(N, double**);
		for(i = 0; i < N; i++) {
			*(jtps + i) = Calloc( T[i], double* );
			for(t = 0; t < T[i]; t++)
				*(*(jtps + i) + t) = Calloc(CRF_params->S * CRF_params->S, double);
		}
	}

	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = T[0];
	for(i = 1; i < N; i++)
		if(maxT < T[i])
			maxT = T[i];

	double **forward_vars = Calloc(maxT, double*);
	double **backward_vars = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(CRF_params->S, double);
		backward_vars[t] = Calloc(CRF_params->S, double);
	}
	
	// pointers to the observed class and covariates at all time points for each subject
	// the observed class
	int *a[N];
	double *X[N];

	for(i = 0; i < N; i++) {
		a[i] = INTEGER(getListElement(VECTOR_ELT(observedCRFsp, i), "y"));
		X[i] = REAL(getListElement(VECTOR_ELT(observedCRFsp, i), "X"));
	}
	
	// DO THE CALCULATION OF THE GRADIENT
	// calculate updated log_obs_probs
	update_log_obs_probs_all_subs_given_state_parametric_CRF(CRF_params, N, T, REAL(betasp), INTEGER(colsp), C,
		X, prev_log_obs_probs_by_state, log_obs_probs_by_state);

	// calculate the marginal probabilities of each class at all t based on all models combined
	for(i = 0; i < N; i++) {
		calc_forward_vars_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), forward_vars);
		calc_backward_vars_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), backward_vars);

		calc_ps_logs(CRF_params, *(T + i), forward_vars, backward_vars, *(ps + i));
		if(diff_wrt_omegas)
			calc_jtps_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), forward_vars, backward_vars, *(jtps + i));
	}

	tempdblptr1 = REAL(retval);

	// partial derivatives wrt omegas
	if(diff_wrt_omegas) {
		if(reduced_log_mat_parameterization) {
			*(tempdblptr1) = partial_CRF_log_lik_wrt_omega_reduced_parameterization(CRF_params, N, T, a, jtps);
		} else {
			for(r = 0; r < CRF_params->S - 1; r++)
				for(s = 0; s < CRF_params->S; s++)
					*(tempdblptr1 + s + r * CRF_params->S) = partial_CRF_log_lik_wrt_omega_r_s(CRF_params, N, T, a, jtps, r, s);

			r = CRF_params->S - 1;
			for(s = 0; s < CRF_params->S - 1; s++)
				*(tempdblptr1 + s + r * CRF_params->S) = partial_CRF_log_lik_wrt_omega_r_s(CRF_params, N, T, a, jtps, r, s);
		}
	}

	// partial derivatives wrt betas
	int beta_diff_start_ind;
	if(diff_wrt_omegas) {
		if(reduced_log_mat_parameterization) {
			beta_diff_start_ind = 1;
		} else {
			beta_diff_start_ind = CRF_params->S * CRF_params->S - 1;
		}
	} else {
		beta_diff_start_ind = 0;
	}

	for(s = 0; s < Sm1; s++) {
		for(d = 0; d < num_betas_diff; d++) {
			*(tempdblptr1 + beta_diff_start_ind + s*num_betas_diff + d) =
				partial_parametric_CRF_log_lik_wrt_beta_s_d_given_obs_probs(CRF_params, N, T, a, X, ps, s,
					*(beta_diff_cols + d) - 1);
		}
	}

	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
		Free(backward_vars[t]);
	}
	Free(forward_vars);
	Free(backward_vars);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(prev_log_obs_probs_by_state + i) + t));
			Free(*(*(log_obs_probs_by_state + i) + t));
			Free(*(*(ps + i) + t));
		}
		Free(*(prev_log_obs_probs_by_state + i));
		Free(*(log_obs_probs_by_state + i));
		Free(*(ps + i));
	}
	Free(prev_log_obs_probs_by_state);
	Free(log_obs_probs_by_state);
	Free(ps);

	if(diff_wrt_omegas) {
		for(i = 0; i < N; i++) {
			for(t = 0; t < T[i]; t++)
				Free(*(*(jtps + i) + t));
			Free(*(jtps + i));
		}
		Free(jtps);
	}


	UNPROTECT(nRprotect);

	return retval;
}



double partial_CRF_log_lik_wrt_omega_r_s(CRF *CRF_params, int N, int *T, int **a, double ***jtps, int r, int s) {
	signedlog *resultsl = new_signedlog(-INFINITY, 1);
	double temp;

	for(int i = 0; i < N; i++)
		for(int t = 1; t < *(T + i); t++)
			if(*(*(a + i) + t - 1) - 1 == r && *(*(a + i) + t) - 1 == s) {
				temp = logspace_sub(0, *(*(*(jtps + i) + t) + r * CRF_params->S + s));
				if(isnan(temp))
					temp = -INFINITY;
				logspace_add_signedlogs(resultsl->slval, resultsl->slsign, temp, 1, &(resultsl->slval), &(resultsl->slsign));
			} else {
				logspace_add_signedlogs(resultsl->slval, resultsl->slsign, *(*(*(jtps + i) + t) + r * CRF_params->S + s), -1, &(resultsl->slval), &(resultsl->slsign));
			}

	temp = signedlog_to_dbl(resultsl);

	free_signedlog(resultsl);

	return(temp);
}



double partial_CRF_log_lik_wrt_omega_reduced_parameterization(CRF *CRF_params, int N, int *T, int **a, double ***jtps) {
	signedlog *resultsl = new_signedlog(-INFINITY, 1);
	double temp;

	for(int i = 0; i < N; i++)
		for(int t = 1; t < *(T + i); t++) {
			temp = -INFINITY;
			for(int s = 0; s < CRF_params->S; s++)
				temp = logspace_add(temp, *(*(*(jtps + i) + t) + s * CRF_params->S + s));

			if(*(*(a + i) + t - 1) == *(*(a + i) + t)) {
				temp = logspace_sub(0, temp);
				if(isnan(temp))
					temp = -INFINITY;
				logspace_add_signedlogs(resultsl->slval, resultsl->slsign, temp, 1, &(resultsl->slval), &(resultsl->slsign));
			} else {
				logspace_add_signedlogs(resultsl->slval, resultsl->slsign, temp, -1, &(resultsl->slval), &(resultsl->slsign));
			}
		}

	temp = signedlog_to_dbl(resultsl);

	free_signedlog(resultsl);

	return(temp);
}



SEXP CRF_gradient_log_lik_wrt_omegas_given_obs_probs_R_interface(SEXP CRF_Xptr, SEXP observedCRFsp, SEXP combined_log_obs_probs_by_statep, SEXP reduced_trans_mat_parameterizationp) {
	int i, t, r, s, maxT, nRprotect = 0, N = length(observedCRFsp), reduced_trans_mat_parameterization = *INTEGER(reduced_trans_mat_parameterizationp);
	int T[N];
	double *tempdblptr1;

	// READ VARIOUS PARAMETER VALUES INTO C DATA STRUCTURES

	// CRF_params contains the parameters for the CRF.
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);

	for(i = 0; i < N; i++)
		T[i] = *(INTEGER(getAttrib(getListElement(VECTOR_ELT(observedCRFsp, i), "X"), R_DimSymbol)));

	// a variable to contain the return value
	SEXP retval;

	if(reduced_trans_mat_parameterization) {
		retval = PROTECT(allocVector(REALSXP, 1));
	} else {
		retval = PROTECT(allocVector(REALSXP, CRF_params->S * CRF_params->S - 1));
	}
	nRprotect++;

	
	double ***log_obs_probs_by_state = Calloc(N, double**);
	double ***jtps = Calloc(N, double**);
	for(i = 0; i < N; i++) {
		*(log_obs_probs_by_state + i) = Calloc( T[i], double* );
		*(jtps + i) = Calloc( T[i], double* );
		for(t = 0; t < T[i]; t++) {
			*(*(log_obs_probs_by_state + i) + t) = Calloc(CRF_params->S, double);
			*(*(jtps + i) + t) = Calloc(CRF_params->S * CRF_params->S, double);
		}
		
		tempdblptr1 = REAL(VECTOR_ELT(combined_log_obs_probs_by_statep, i));
		for(t = 0; t < T[i]; t++) {
			for(s = 0; s < CRF_params->S; s++) {
				*(*(*(log_obs_probs_by_state + i) + t) + s) = *(tempdblptr1 + t + s * T[i]);
			}
		}
	}
	
	// create space to store forward and backward variables, observation probabilities, and combined class probabilities.
	// it is convenient for the forward vars to be stored in a row vector, backward vars in a col vector,
	// and probabilities given states in a col vector
	// We allocate enough space to store the variables for the subject with the longest observation sequence
	// and re-use that space for all subjects
	maxT = T[0];
	for(i = 1; i < N; i++)
		if(maxT < T[i])
			maxT = T[i];

	double **forward_vars = Calloc(maxT, double*);
	double **backward_vars = Calloc(maxT, double*);
	for(t = 0; t < maxT; t++) {
		forward_vars[t] = Calloc(CRF_params->S, double);
		backward_vars[t] = Calloc(CRF_params->S, double);
	}
	
	// pointers to the observed class at all time points for each subject
	int *a[N];

	for(i = 0; i < N; i++)
		a[i] = INTEGER(getListElement(VECTOR_ELT(observedCRFsp, i), "y"));
	
	// DO THE CALCULATION OF THE GRADIENT
	// calculate the probabilities of each class at all t based on all models combined
	for(i = 0; i < N; i++) {
		calc_forward_vars_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), forward_vars);
		calc_backward_vars_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), backward_vars);

		calc_jtps_logs(CRF_params, *(T + i), *(log_obs_probs_by_state + i), forward_vars, backward_vars, *(jtps + i));
	}

	tempdblptr1 = REAL(retval);
	// gradient wrt omegas
	if(reduced_trans_mat_parameterization) {
		*(tempdblptr1) = partial_CRF_log_lik_wrt_omega_reduced_parameterization(CRF_params, N, T, a, jtps);
	} else {
		for(r = 0; r < CRF_params->S - 1; r++)
			for(s = 0; s < CRF_params->S; s++)
				*(tempdblptr1 + s + r * CRF_params->S) = partial_CRF_log_lik_wrt_omega_r_s(CRF_params, N, T, a, jtps, r, s);

		r = CRF_params->S - 1;
		for(s = 0; s < CRF_params->S - 1; s++)
			*(tempdblptr1 + s + r * CRF_params->S) = partial_CRF_log_lik_wrt_omega_r_s(CRF_params, N, T, a, jtps, r, s);
	}

	// FREE MEMORY
	for(t = 0; t < maxT; t++) {
		Free(forward_vars[t]);
		Free(backward_vars[t]);
	}
	Free(forward_vars);
	Free(backward_vars);

	for(i = 0; i < N; i++) {
		for(t = 0; t < T[i]; t++) {
			Free(*(*(log_obs_probs_by_state + i) + t));
			Free(*(*(jtps + i) + t));
		}
		Free(*(log_obs_probs_by_state + i));
		Free(*(jtps + i));
	}
	Free(log_obs_probs_by_state);
	Free(jtps);

	UNPROTECT(nRprotect);

	return retval;
}





////////////////////////////////////////////////////////////////////////////////
// Make C functions available in R namespace
////////////////////////////////////////////////////////////////////////////////

R_CallMethodDef callMethods[] =
{
    {"initialize_new_CRF_C", (DL_FUNC)&initialize_new_CRF, 4},
    {"set_CRF_pi_from_log_C", (DL_FUNC)&set_CRF_pi_from_log, 2},
    {"set_CRF_trans_matrix_from_log_C", (DL_FUNC)&set_CRF_trans_matrix_from_log, 2},
	{"parametric_CRF_update_log_lik_C", (DL_FUNC)&parametric_CRF_update_log_lik_R_interface, 6},
	{"CRF_log_lik_given_obs_probs_C", (DL_FUNC)&CRF_log_lik_given_obs_probs_R_interface, 3},
	{"parametric_CRF_update_gradient_log_lik_C", (DL_FUNC)&parametric_CRF_update_gradient_log_lik_R_interface, 8},
	{"CRF_gradient_log_lik_wrt_omegas_given_obs_probs_C", (DL_FUNC)&CRF_gradient_log_lik_wrt_omegas_given_obs_probs_R_interface, 4},
	{"calc_log_McShane_class_probs_given_log_static_class_probs_C", (DL_FUNC)&calc_log_McShane_class_probs_given_log_static_class_probs_R_interface, 5},
	{"marginal_CRF_log_lik_given_obs_probs_C", (DL_FUNC)&marginal_CRF_log_lik_given_obs_probs_R_interface, 4},
	{"get_dbl_max_C", (DL_FUNC)&get_dbl_max, 0},
	{"logspace_add_C", (DL_FUNC)&logspace_add_R_C_interface, 2},
	{"logspace_sub_C", (DL_FUNC)&logspace_sub_R_C_interface, 2},
    {"logspace_sum_matrix_rows_C", (DL_FUNC)&logspace_sum_matrix_rows, 3},
    {"logspace_sub_matrix_rows_C", (DL_FUNC)&logspace_sub_matrix_rows, 2},
    {"dMVN_multiobs_C", (DL_FUNC)&dMVN_multiobs_R_interface, 7},
    {"dMVN_Prec_multiobs_C", (DL_FUNC)&dMVN_Prec_multiobs_R_interface, 7},
    {"dMVN_DiagPrec_multiobs_C", (DL_FUNC)&dMVN_DiagPrec_multiobs_R_interface, 7},
    {"dGMM_sameDiagPrec_C", (DL_FUNC)&dGMM_sameDiagPrec_R_Interface, 8},
    {"dGMM_DiagPrec_C", (DL_FUNC)&dGMM_DiagPrec_R_Interface, 8},
    {"dGMM_samePrec_C", (DL_FUNC)&dGMM_samePrec_R_Interface, 8},
    {"dGMM_Prec_C", (DL_FUNC)&dGMM_Prec_R_Interface, 8},
    {"dGMM_sameSigma_C", (DL_FUNC)&dGMM_sameSigma_R_Interface, 8},
    {"dGMM_C", (DL_FUNC)&dGMM_R_Interface, 8},
	{NULL,NULL, 0}
};

void R_init_PACwithDDM(DllInfo *dll)
{
	R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}
