/*
 ============================================================================
 Name        : utility.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : utility functions for numeric calculations and interacting with R
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
#include <float.h>

#include "utility.h"


SEXP get_dbl_max() {
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	*(REAL(retval)) = DBL_MAX;
	
	UNPROTECT(1);
	return retval;
}

SEXP getListElement(SEXP list, const char *str) {
	// Return the element with name matching str from the list provided.
	// Based on http://cran.r-project.org/doc/manuals/R-exts.html#Handling-lists
	//
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	for (R_len_t i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}

	return elmt;
}


SEXP logspace_add_R_C_interface(SEXP log_x, SEXP log_y) {
	// An interface to R's C function logspace_add
	// Computes log(exp(log_x) - exp(log_y))
	// Returns -INFINITY instead of NaN if both log_x and log_y are -INFINITY
	
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));
	
	*(REAL(retval)) = logspace_add_safe(lx, ly);
	
	UNPROTECT(1);
	return retval;
}


double logspace_add_safe(double log_x, double log_y) {
	// An interface to R's C function logspace_add that
	// Computes log(exp(log_x) + exp(log_y))
	// Returns -INFINITY instead of NaN if both log_x and log_y are -INFINITY
	
	if(log_x == R_NegInf && log_y == R_NegInf) {
		return(R_NegInf);
	} else {
		return(logspace_add(log_x, log_y));
	}
}


SEXP logspace_sub_R_C_interface(SEXP log_x, SEXP log_y) {
	// An interface to R's C function logspace_sub
	// Computes log(exp(log_x) - exp(log_y))
	
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));

	*(REAL(retval)) = logspace_sub(lx, ly);
	
	UNPROTECT(1);
	return retval;
}

SEXP logspace_sum_matrix_rows(SEXP Xp, SEXP N_rowp, SEXP N_colp) {
	int i, j, n_row = *INTEGER(N_rowp), n_col = *INTEGER(N_colp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = *(X + i);
	}

	for(j = 1; j < n_col; j++) {
		for(i = 0; i < n_row; i++) {
			if(!(*(dblptr + i) == R_NegInf && *(X + i + j*n_row) == R_NegInf))
				*(dblptr + i) = logspace_add(*(dblptr + i), *(X + i + j*n_row));
		}
	}

	UNPROTECT(1);
	return retval;
}

SEXP logspace_sub_matrix_rows(SEXP Xp, SEXP N_rowp) {
	int i, n_row = *INTEGER(N_rowp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = logspace_sub(*(X + i), *(X + i + n_row));
	}

	UNPROTECT(1);
	return retval;
}


int intmin(int a, int b) {
	if(a < b) {
		return a;
	} else {
		return b;
	}
}


int intmax(int a, int b) {
	if(a > b) {
		return a;
	} else {
		return b;
	}
}


double calc_matrix_determinant(double* matrixp, int length, int retlog) {
	/* INPUTS:
	 * matrixp is a pointer to the first entry of the matrix
	 * length is the number of rows (and columns) in the matrix
	 *
	 * OUTPUTS:
	 * Returns the determinant of the matrix
	 */
	int i;

	double mat[length*length];

	double logdet = 0;

	int eigres, worklen = 4*length;
	double work[worklen], evals[length];
    char jobz[1], uplo[1];
	jobz[0] = 'N';
	uplo[0] = 'L';

	for(i = 0; i < length*length; i++) {
		mat[i] = *(matrixp + i);
	}

	// dsyev computes the eigenvalues of a symmetric matrix
	F77_CALL(dsyev)(jobz, uplo, &length, (double*)&mat, &length, (double*)&evals, (double*)&work, &worklen, &eigres);

	// the logarithm of the determinant is the sum of log-eigenvalues
	for(i = 0; i < length; i ++) {
		logdet += log(evals[i]);
	}

	if(retlog) {
		return logdet;
	} else {
		return exp(logdet);
	}
}


void matrix_multiply(double *A, double *B, double *result, int M, int K, int N) {
	// compute the matrix product A B and store in result
	// A is an M by K matrix, B is K by N matrix
	char transa[1], transb[1];
	transa[0] = 'N', transb[0] = 'N';
	double alpha = 1, beta = 0;

	if(result == A || result == B) {
		double *temp = Calloc(M * N, double);
		F77_CALL(dgemm)(transa, transb, &M, &N, &K, &alpha, A, &M, B, &K, &beta, temp, &M);
		for(int i = 0; i < M * N; i++) {
			*(result + i) = *(temp + i);
		}
		Free(temp);
	} else {
		F77_CALL(dgemm)(transa, transb, &M, &N, &K, &alpha, A, &M, B, &K, &beta, result, &M);
	}
}


void transform_p_to_f(double *ps, double *fs, int len) {
	int i, j, info, adjlen = len-1, intone = 1;
	double tempmat[adjlen*adjlen];
	//	double temp;
	int junk[adjlen];

	// adjust ps so that they are between 10^-8 and 1 - len*10^-8, to eliminate numerical problems
	// this works ok for a small number of categories.  do something more real for a larger number of categories.
   	for(i = 0; i < adjlen; i++) {
		if(*(ps + i) > 1 - len*0.00000001) {
			*(ps + i) = 1 - len*0.00000001;
		} else if(*(ps + i) < 0.00000001) {
			*(ps + i) = 0.00000001;
		}
	}

	// be sure that the ps sum to 1
	//	temp = 0.0;
	//	for(i = 0; i < len; i++)
	//	  temp += *(ps + i);

	//	for(i = 0; i < len; i++)
	//	  *(ps + i) /= temp;

	// tempmat = I_len-1 - P (1_len-1)'
	for(i = 0; i < adjlen; i++)
		for(j = 0; j < adjlen; j++)
			tempmat[i + j*adjlen] = -1 * *(ps + i);

	for(i = 0; i < adjlen; i++)
		tempmat[i + i*adjlen] += 1;

	// store ps in f: the argument to the matrix multiply function gets modified to contain the result
	for(i = 0; i < adjlen; i++)
		*(fs + i) = *(ps + i);

	// solve expfs = tempmat^-1 ps
	F77_CALL(dgesv)(&adjlen, &intone, (double*)&tempmat, &adjlen, (int*)&junk, fs, &adjlen, &info);

	// compute fs from expfs
	for(i = 0; i < adjlen; i++)
		*(fs + i) = log(*(fs + i));

	// the f function for the final category is 0
	*(fs + len - 1) = 0;
}
