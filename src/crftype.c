/*
 ============================================================================
 Name        : crftype.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : Type definitions and functions for CRFs
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <math.h>

#include "crftype.h"
#include "utility.h"


SEXP initialize_new_CRF(SEXP Dp, SEXP Sp, SEXP log_pi, SEXP log_trans_mat) {
	// Dp is the number of covariates observed at each time point
	// Sp is the number of states in the state space
	// log_pi is log(pi), where pi is proportional to the probabilities of the initial state distribution
	// log_trans_mat is log(trans_mat), where trans_mat is proportional to the transition probabilities

	CRF *CRF_params = Calloc(1, CRF);

	CRF_params->D = *(INTEGER(Dp));
	CRF_params->S = *(INTEGER(Sp));

	CRF_params->log_pi = Calloc( CRF_params->S, double);
	CRF_params->log_trans_matrix = Calloc( CRF_params->S * CRF_params->S, double);
	
	SEXP CRF_Xptr = CRF_RXptr(CRF_params);
	
	set_CRF_pi_from_log(CRF_Xptr, log_pi);
	set_CRF_trans_matrix_from_log(CRF_Xptr, log_trans_mat);
	
	return(CRF_Xptr);
}


SEXP set_CRF_pi_from_log(SEXP CRF_Xptr, SEXP log_pi) {
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);
	
	double *tempdbl = REAL(log_pi);
	for(int i = 0; i < CRF_params->S; i++)
		*(CRF_params->log_pi + i) = *(tempdbl + i);

	return(R_NilValue);
}

SEXP set_CRF_trans_matrix_from_log(SEXP CRF_Xptr, SEXP log_trans_mat) {
	CRF *CRF_params = CRF_from_RXptr(CRF_Xptr);

	double *tempdbl = REAL(log_trans_mat);
	for(int i = 0; i < CRF_params->S * CRF_params->S; i++)
		*(CRF_params->log_trans_matrix + i) = *(tempdbl + i);
	
	return(R_NilValue);
}

void deconstruct_CRF(CRF *CRF_params) {
	Free(CRF_params->log_pi);
	Free(CRF_params->log_trans_matrix);
}

void _finalize_CRF_RXptr(SEXP ext) {
	if(!R_ExternalPtrAddr(ext))
		return;
	
	CRF *CRF_params = CRF_from_RXptr(ext);
	deconstruct_CRF(CRF_params);
	Free(CRF_params);
	
	R_ClearExternalPtr(ext);
}

SEXP CRF_RXptr(CRF *CRF_params) {
	SEXP ext = PROTECT(R_MakeExternalPtr(CRF_params, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext, _finalize_CRF_RXptr, TRUE);
	UNPROTECT(1);
	
	return ext;
}

CRF* CRF_from_RXptr(SEXP CRF_Xptr) {
	return (CRF *)R_ExternalPtrAddr(CRF_Xptr);
}
