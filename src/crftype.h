#ifndef CRFTYPE_H
#define CRFTYPE_H

#include <Rdefines.h>

/*
 ============================================================================
 Name        : crftype.h
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : Type definitions and functions for CRFs
 ============================================================================
 */

typedef struct {
	int D; // Number of covariates
	int S; // State space size
	double *log_pi;
	double *log_trans_matrix;
} CRF;


SEXP initialize_new_CRF(SEXP Dp, SEXP Sp, SEXP log_pi, SEXP log_trans_mat);

SEXP set_CRF_pi_from_log(SEXP CRF_ptr, SEXP log_pi);

SEXP set_CRF_trans_matrix_from_log(SEXP CRF_Xptr, SEXP log_trans_mat);

void deconstruct_CRF(CRF *CRF_params);

void _finalize_CRF_RXptr(SEXP ext);

SEXP CRF_RXptr(CRF *CRF_params);

CRF* CRF_from_RXptr(SEXP CRF_Xptr);

#endif
