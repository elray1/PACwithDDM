/*
 ============================================================================
 Name        : scinotFunctions.h
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : Type definition and functions for scientific notation matrices
 ============================================================================
 */

typedef struct {
	double *coef;
	int *expon;
	int nrow;
	int ncol;
} scinot_matrix;

scinot_matrix alloc_scinot(int nrow, int ncol);

scinot_matrix duplicate_scinot(scinot_matrix *mat);

void free_scinot(scinot_matrix *mat);
	
void Rprintf_scinot(scinot_matrix *mat);

void to_scinot(double *x, int length, double *coef, int *expon);

void log_to_scinot(double *x, int length, double *coef, int *expon);

void scinot_to_double(scinot_matrix *source, double* target);

void copy_scinot(scinot_matrix *source, scinot_matrix *target);

void insert_scinot(scinot_matrix *mat1, scinot_matrix *mat2, int rowind, int colind);

void normalize_scinot(scinot_matrix *mat, scinot_matrix *target);

void center_scinot(scinot_matrix *mat, scinot_matrix *target);

void adjust_expon_scinot(int length, double *coef, int *expon);

void sum_all_entries_scinot(scinot_matrix *mat, scinot_matrix *target);

void sum_cols_scinot(scinot_matrix *mat, scinot_matrix *target);

void sum_rows_scinot(scinot_matrix *mat, scinot_matrix *target);

void reciprocal_scinot(scinot_matrix *mat);
	
void add_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target);

void add_one_to_one_scinot(scinot_matrix *mat1, int rowind1, int colind1, scinot_matrix *mat2, int rowind2, int colind2, scinot_matrix *target);

void add_to_one_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target, int rowind, int colind);

void add_to_all_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target);

void mult_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target);

void mult_one_to_one_scinot(scinot_matrix *mat1, int rowind1, int colind1, scinot_matrix *mat2, int rowind2, int colind2, scinot_matrix *target);
	
void mult_to_one_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target, int rowind, int colind);

void mult_to_all_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target);

void matrix_mult_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target);

void log_sum_all_entries_exp_scinot(scinot_matrix *source, scinot_matrix *target);

void log_sum_all_entries_exp_scinotDBG(scinot_matrix *source, scinot_matrix *target);

void exp_scinot(scinot_matrix *source, scinot_matrix *target);

void exp_scinotDBG(scinot_matrix *source, scinot_matrix *target);

void log_scinot(scinot_matrix *source, scinot_matrix *target);

void sqrt_scinot(scinot_matrix *source, scinot_matrix *target);
