/* 
 ============================================================================
 Name        : utilityFunctions.h
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : utility functions for numeric calculations and interacting with R
 ============================================================================
 */

#include <Rdefines.h>

SEXP get_dbl_max();

SEXP getListElement(SEXP list, const char *str);

SEXP logspace_add_R_C_interface(SEXP x, SEXP y);

double logspace_add_safe(double log_x, double log_y);

SEXP logspace_sub_R_C_interface(SEXP log_x, SEXP log_y);

SEXP logspace_sum_matrix_rows(SEXP Xp, SEXP N_rowp, SEXP N_colp);

SEXP logspace_sub_matrix_rows(SEXP Xp, SEXP N_rowp);

int intmin(int a, int b);

int intmax(int a, int b);

double calc_matrix_determinant(double* matrixp, int length, int retlog);

void matrix_multiply(double *A, double *B, double *result, int M, int K, int N);

void transform_p_to_f(double *ps, double *fs, int len);
