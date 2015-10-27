/*
 ============================================================================
 Name        : scinotFunctions.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : Hello World in C, Ansi-style
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

#include <math.h>

#include "scinot.h"
#include "utility.h"


scinot_matrix alloc_scinot(int nrow, int ncol) {
	scinot_matrix retval;
	retval.nrow = nrow;
	retval.ncol = ncol;
	retval.coef = Calloc(nrow*ncol, double);
	retval.expon = Calloc(nrow*ncol, int);
//	retval.coef = (double*)calloc(nrow*ncol, sizeof(double));
//	retval.expon = (int*)calloc(nrow*ncol, sizeof(int));
	return retval;
}
/*scinot_matrix alloc_scinot(int nrow, int ncol) {
	scinot_matrix retval;
	retval.nrow = nrow;
	retval.ncol = ncol;
	retval.coef = (double*)R_alloc(nrow*ncol, sizeof(double));
	retval.expon = (int*)R_alloc(nrow*ncol, sizeof(int));
	return retval;
}*/

scinot_matrix duplicate_scinot(scinot_matrix *mat) {
	int i;
	scinot_matrix retval = alloc_scinot(mat->nrow, mat->ncol);

	for(i = 0; i < mat->nrow * mat->ncol; i++) {
		*(retval.coef + i) = *(mat->coef + i);
		*(retval.expon + i) = *(mat->expon + i);
	}

	return retval;
}


void free_scinot(scinot_matrix *mat) {
	Free(mat->coef);
	Free(mat->expon);
//	free(mat->coef);
//	free(mat->expon);
}

void Rprintf_scinot(scinot_matrix *mat) {
	int i, j;

	for(i = 0; i < mat->nrow; i++) {
		Rprintf("\n");
		for(j = 0; j < mat->ncol; j++) {
			Rprintf("%lf * 10^%i, ", *(mat->coef + i + j*mat->nrow), *(mat->expon + i + j*mat->nrow));
		}
	}
}

void to_scinot(double *x, int length, double *coef, int *expon) {
/* calculates x in sci notation with a separate power for each entry so that
 * (10^(-expon_ij))*x_ij is between +/-0.1 and +/-1.
 * If x_ij == 0, the exponent is set to be 0.
 * expon is 10s power, coef is coefficient
 */
	int i;
	double temp = 0.0;

	for(i = 0; i < length; i++) {
		temp = fabs(*(x + i));

		if(temp > DBL_EPSILON) {
			// find appropriate exponent based on the magnitude of x_ij
			*(expon + i) = ceil(log10(temp));
		} else {
			// if x_ij is (approximately) 0, cancel out the exponent
			*(expon + i) = 0;
		}

		// adjust coefficients to match the exponent
		*(coef + i) = *(x + i) * pow(10, - *(expon + i));
	}
}


void log_to_scinot(double *logx, int length, double *coef, int *expon) {
/* calculates exp(logx) in sci notation at a power so that
 * (10^(-expon))*sum(exp(logx)) is between +/-0.1 and +/-1.
 * Note that each element of the result is given the same order.
 * expon is 10s power, coef is coefficient
 */
	int i;
	double temp = 0.0;
	double log10 = log(10);

	// convert each element of logx individually
	for(i = 0; i < length; i++) {
		if(*logx == -INFINITY) {
		// the exponent for x_i
			*(expon + i) = 0;
//		Rprintf("\nexpon: %i", *(expon + i));
		// the coefficient for x_i
			*(coef + i) = 0;
		} else {
		// temp = log10(x_i)
//		Rprintf("\nlog(x): %lf", *(logx + i));
			temp = *(logx + i)/log10;
//		Rprintf("\nlog_10(x): %lf", temp);
		// the exponent for x_i
			*(expon + i) = ceil(temp);
//		Rprintf("\nexpon: %i", *(expon + i));
		// the coefficient for x_i
			*(coef + i) = pow(10, temp - *(expon + i));
		}
	}
}


void scinot_to_double(scinot_matrix *source, double* target) {
/* convert source to a double, store in target.
 * it is assumed that target has the same number of elements as source
 */
	int i;

	for(i = 0; i < source->nrow * source->ncol; i++) {
		*(target + i) = *(source->coef + i) * pow(10, *(source->expon + i));
	}
}


void copy_scinot(scinot_matrix *source, scinot_matrix *target) {
/* copy source to target
 * it is assumed that source and target have the same dimensions.
 */
	int i;

	for(i = 0; i < source->nrow * source->ncol; i++) {
		*(target->coef + i) = *(source->coef + i);
		*(target->expon + i) = *(source->expon + i);
	}
}


void insert_scinot(scinot_matrix *mat1, scinot_matrix *mat2, int rowind, int colind) {
/* mat1 is a matrix in scientific notation, mat2 is a single value in scientific notation
 * insert mat2 into mat1 at the indices indicated by rowind, colind
 * it is assumed that mat2 is a 1 by 1 matrix, and rowind < mat1->nrow and colind < mat1->ncol
 */

	*(mat1->coef + rowind + colind * mat1->nrow) = *(mat2->coef);
	*(mat1->expon + rowind + colind * mat1->nrow) = *(mat2->expon);
}


void adjust_expon_scinot(int length, double *coef, int *expon) {
/* adjust the exponent so that each coef is between +/-0.1 and +/-1
 */
	int i, exponadj;
	double temp;

	for(i = 0; i < length; i++) {
		temp = fabs(*(coef + i));
	
		if(temp > DBL_EPSILON) {
			// if the coefficient is not (approximately) 0, consider making an adjustment to the exponent.
			exponadj = ceil(log10(temp));

			// adjust the exponent and coefficients to match the exponent
			// for example, if exponadj == -3, the coefficients are too small.
			// multiply all coefficients by 10^3 to make them larger and
			// subtract 3 from the exponent to counteract.
			if(exponadj != 0) {
				*(expon + i) += exponadj;
				*(coef + i) = *(coef + i) * pow(10, -exponadj);
			}
		} else {
			*(expon + i) = 0;
		}
	}
}


void normalize_scinot(scinot_matrix *mat, scinot_matrix *target) {
/* normalize mat so that the sum of all entries is 1
 * Basically, divide the matrix by sum(coef*10^expon) = sum(coef)*10^expon
 */
	scinot_matrix temp = alloc_scinot(1, 1);

	sum_all_entries_scinot(mat, &temp);
//	Rprintf("sum of all entries is %lf * 10^%i", *(temp.coef), *(temp.expon));

	// it's only possible to scale the matrix so that entries sum to 1 if the entries don't sum to 0.
	if(*(temp.coef) != 0.0) {
		reciprocal_scinot(&temp);
//		Rprintf("inverse sum of all entries is %lf * 10^%i", *(temp.coef), *(temp.expon));
		mult_to_all_scinot(mat, &temp, target);
//		Rprintf("normalized matrix is %lf * 10^%i", *(temp.coef), *(temp.expon));
	}

	free_scinot(&temp);
}


void center_scinot(scinot_matrix *mat, scinot_matrix *target) {
/* center mat by subtracting the mean from every value
 */
	scinot_matrix temp1 = alloc_scinot(1, 1), temp2 = alloc_scinot(1, 1);

	// -1 * mean value
	sum_all_entries_scinot(mat, &temp1);
//Rprintf("In centering: sum all values is %lf * 10^%i", *(temp1.coef), *(temp1.expon));
	*(temp2.coef) = -1 * mat->nrow * mat->ncol;
	*(temp2.expon) = 0;
	reciprocal_scinot(&temp2);

	mult_scinot(&temp1, &temp2, &temp1);
//Rprintf("\nIn centering: mean value is %lf * 10^%i", *(temp1.coef), *(temp1.expon));

	// subtract mean from every entry of mat
	add_to_all_scinot(mat, &temp1, target);
//Rprintf("\nExit center");
	free_scinot(&temp1);
	free_scinot(&temp2);
}


void sum_all_entries_scinot(scinot_matrix *mat, scinot_matrix *target) {
/* target gets the sum of all entries of mat
 * it is assumed that target has 1 row 1 column
 */
	int i, newexpon, minexpon;

	// find the smallest exponent in the matrix.
	// We do this as a first step in finding the largest exponent for a non-zero number
	// could save time by making this a while loop
	minexpon = *(mat->expon);	
	for(i = 1; i < mat->nrow * mat->ncol; i++) {
		if(*(mat->expon + i) < minexpon)
			minexpon = *(mat->expon + i);
	}

	// find the largest exponent for a non-zero number
	newexpon = minexpon;
	for(i = 0; i < mat->nrow * mat->ncol; i++) {
		if(*(mat->expon + i) > newexpon && fabs(*(mat->coef + i)) > DBL_EPSILON)
			newexpon = *(mat->expon + i);
	}

	// sum, adjusting each term as appropriate to get the right final exponent
	*(target->coef) = 0;
	*(target->expon) = newexpon;
	for(i = 0; i < mat->nrow * mat->ncol; i++) {
		if(fabs(*(mat->coef + i)) > DBL_EPSILON)
			*(target->coef) = *(target->coef) + *(mat->coef + i) * pow(10, *(mat->expon + i) - newexpon);
	}

	// although newexpon is the approximately correct exponent,
	// the summation of (scaled) coefficients may have thrown things off
	adjust_expon_scinot(1, target->coef, target->expon);
}


void sum_cols_scinot(scinot_matrix *mat, scinot_matrix *target) {
/* target gets the column sums of mat
 * it is assumed that target has 1 row and the same number of columns as mat
 */
	int i, j, newexpon, minexpon;

	// find the smallest exponent in the matrix.
	// We do this as a first step in finding the largest exponent for a non-zero number within each column
	// could save time by making this a while loop
	minexpon = *(mat->expon);	
	for(i = 1; i < mat->nrow * mat->ncol; i++) {
		if(*(mat->expon + i) < minexpon)
			minexpon = *(mat->expon + i);
	}

	// for each column, ...
	for(j = 0; j < mat->ncol; j++) {
		// find the largest exponent for a non-zero number within that column
		newexpon = minexpon;
		for(i = 0; i < mat->nrow; i++) {
			if(*(mat->expon + i + j*mat->nrow) > newexpon && fabs(*(mat->coef + i + j*mat->nrow)) > DBL_EPSILON)
				newexpon = *(mat->expon + i + j*mat->nrow);
		}

		// sum the terms in the column, adjusting each term as appropriate to get the right final exponent
		*(target->coef + j) = 0;
		*(target->expon + j) = newexpon;
		for(i = 0; i < mat->nrow; i++) {
			if(fabs(*(mat->coef + i + j*mat->nrow)) > DBL_EPSILON)
				*(target->coef + j) = *(target->coef + j) + *(mat->coef + i + j*mat->nrow) * pow(10, *(mat->expon + i + j*mat->nrow) - newexpon);
		}
	}

	// although newexpon is the approximately correct exponent,
	// the summation of (scaled) coefficients may have thrown things off
	adjust_expon_scinot(1 * mat->ncol, target->coef, target->expon);
}


void sum_rows_scinot(scinot_matrix *mat, scinot_matrix *target) {
/* target gets the row sums of mat
 * it is assumed that target has 1 row and the same number of rows as mat
 */
	int i, j, newexpon, minexpon;

	// find the smallest exponent in the matrix.
	// We do this as a first step in finding the largest exponent for a non-zero number within each column
	// could save time by making this a while loop
	minexpon = *(mat->expon);	
	for(i = 1; i < mat->nrow * mat->ncol; i++) {
		if(*(mat->expon + i) < minexpon)
			minexpon = *(mat->expon + i);
	}

	// for each row, ...
	for(i = 0; i < mat->nrow; i++) {
		// find the largest exponent for a non-zero number within that row
		newexpon = minexpon;
		for(j = 0; j < mat->ncol; j++) {
			if(*(mat->expon + i + j*mat->nrow) > newexpon && fabs(*(mat->coef + i + j*mat->nrow)) > DBL_EPSILON)
				newexpon = *(mat->expon + i + j*mat->nrow);
		}

		// sum the terms in the row, adjusting each term as appropriate to get the right final exponent
		*(target->coef + i) = 0;
		*(target->expon + i) = newexpon;
		for(j = 0; j < mat->ncol; j++) {
			if(fabs(*(mat->coef + i + j*mat->nrow)) > DBL_EPSILON)
				*(target->coef + i) = *(target->coef + i) + *(mat->coef + i + j*mat->nrow) * pow(10, *(mat->expon + i + j*mat->nrow) - newexpon);
		}
	}

	// although newexpon is the approximately correct exponent,
	// the summation of (scaled) coefficients may have thrown things off
	adjust_expon_scinot(1 * mat->nrow, target->coef, target->expon);
}


void reciprocal_scinot(scinot_matrix *mat) {
/* invert each element of mat
 */
	int i;

	// 1/(a*10^b) = (1/a) * 10^-b
	for(i = 0; i < mat->nrow * mat->ncol; i++) {
		*(mat->coef + i) = 1 / *(mat->coef + i);
		*(mat->expon + i) *= -1;
	}

	// the inversion of coefficients may have thrown things off
	// (i.e. if coefficients were small, they will now be large)
	adjust_expon_scinot(mat->nrow * mat->ncol, mat->coef, mat->expon);
}


void add_one_term_scinot(double *coef1, int *expon1, double *coef2, int *expon2, double *coeftarget, int *expontarget) {
// add one term.  target may be the same as 1 or 2
// the calculations are based on the fact that if M = max(expon1, expon2),
// coef1*10^expon1 + coef2*10^expon2 = (coef1*10^(expon1-M) + coef2*10^(expon2-M))*10^M
// and M is approximately the right exponent to use since this calculation makes one of the sets of coefficients smaller

	int expadj1, expadj2, newexpon;
	double abscoef1 = fabs(*(coef1)), abscoef2 = fabs(*(coef2));

//	Rprintf("dbl eps is %lf  ", DBL_EPSILON);

	if((abscoef1 > DBL_EPSILON) && (abscoef2 > DBL_EPSILON)) {
		// if both terms are non-zero, add as usual
		newexpon = intmax( *(expon1), *(expon2) );
		expadj1 = *(expon1) - newexpon;
		expadj2 = *(expon2) - newexpon;

		*(coeftarget) = *(coef1) * pow(10, expadj1) + *(coef2) * pow(10, expadj2);
		*(expontarget) = newexpon;
//Rprintf("** %lf, %lf **", abscoef1, abscoef2);
		// although we found the approximately correct exponent,
		// the summation of (scaled) coefficients may have thrown things off
		adjust_expon_scinot(1, coeftarget, expontarget);
	} else if((abscoef1 > DBL_EPSILON) && (abscoef2 <= DBL_EPSILON)) {
		// if the term from mat2 is zero, set equal to the term from mat1
		*(coeftarget) = *(coef1);
		*(expontarget) = *(expon1);
//Rprintf("*1 %lf, %lf 1*", abscoef1, abscoef2);
	} else if((abscoef1 <= DBL_EPSILON) && (abscoef2 > DBL_EPSILON)) {
		// if the term from mat1 is zero, set equal to the term from mat2
		*(coeftarget) = *(coef2);
		*(expontarget) = *(expon2);
//Rprintf("*2 %lf, %lf 2*", abscoef1, abscoef2);
	} else {
		// if the terms from both matrices are zero, set equal to zero
		*(coeftarget) = 0;
		*(expontarget) = 0;
//Rprintf("00 %lf, %lf 00", abscoef1, abscoef2);
	}
}

void add_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target) {
/* add mat1 and mat2, store in target
 * it is assumed that mat1, mat2, and target have the same dimensions.
 */
	int i;

	for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
		add_one_term_scinot(mat1->coef + i, mat1->expon + i, mat2->coef + i, mat2->expon + i, target->coef + i, target->expon + i);
	}
}


void add_one_to_one_scinot(scinot_matrix *mat1, int rowind1, int colind1, scinot_matrix *mat2, int rowind2, int colind2, scinot_matrix *target) {
/* add the specified index of mat2 to the specified index of mat1
 * it is assumed that rowind < mat1->nrow, colind < mat1->ncol, and mat2 is a 1 by 1 matrix
 */
	int i;

	// if mat1 and target don't point to the same object, copy all elements of mat1 to target
	if(mat1 != target) {
		for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
			*(target->coef + i) = *(mat1->coef + i);
			*(target->expon + i) = *(mat1->expon + i);
		}
	}

	add_one_term_scinot(mat1->coef + rowind1 + colind1*mat1->nrow, mat1->expon + rowind1 + colind1*mat1->nrow, mat2->coef + rowind2 + colind2*mat2->nrow, mat2->expon + rowind2 + colind2*mat2->nrow, target->coef + rowind1 + colind1*mat1->nrow, target->expon + rowind1 + colind1*mat1->nrow);
}


void add_to_one_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target, int rowind, int colind) {
/* add mat2 to the specified index of mat1
 * it is assumed that rowind < mat1->nrow, colind < mat1->ncol, and mat2 is a 1 by 1 matrix
 */
	int i;

	// if mat1 and target don't point to the same object, copy all elements of mat1 to target
	if(mat1 != target) {
//		Rprintf("Copying object.\n");
		for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
			*(target->coef + i) = *(mat1->coef + i);
			*(target->expon + i) = *(mat1->expon + i);
		}
	}

	add_one_term_scinot(mat1->coef + rowind + colind*mat1->nrow, mat1->expon + rowind + colind*mat1->nrow, mat2->coef, mat2->expon, target->coef + rowind + colind*mat1->nrow, target->expon + rowind + colind*mat1->nrow);
}


void add_to_all_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target) {
/* add mat2 to every entry of mat1 and store in target
 * it is assumed that mat2 is a 1 by 1 matrix
 */
	int i;

	// the calculations are based on the fact that if M = max(expon1, expon2),
	// coef1*10^expon1 + coef2*10^expon2 = (coef1*10^(expon1-M) + coef2*10^(expon2-M))*10^M
	// and M is approximately the right exponent to use since this calculation makes one of the sets of coefficients smaller
	for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
		add_one_term_scinot(mat1->coef + i, mat1->expon + i, mat2->coef, mat2->expon, target->coef + i, target->expon + i);
	}
}


void mult_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target) {
/* element-wise multiply mat1 and mat2, store in target
 * it is assumed that mat1, mat2, and target have the same dimensions.
 */
	int i;

	// to multiply numbers in scientific notation, we multiply the coefficients and add the exponents.
	for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
		*(target->expon + i) = *(mat1->expon + i) + *(mat2->expon + i);
		*(target->coef + i) = *(mat1->coef + i) * *(mat2->coef + i);
	}

	// although the sum of the exponents is the approximately correct exponent,
	// the multiplication of coefficients may have thrown things off
	adjust_expon_scinot(target->nrow * target->ncol, target->coef, target->expon);
}


void mult_to_one_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target, int rowind, int colind) {
/* multiply element (rowind, colind) of mat1 by mat2
 * it is assumed that mat1 mat2 is 1 by 1, target has the same dimensions as mat1, and
 * rowind < mat1->nrow and colind < mat1->ncol
 */
	int i;

	// if mat1 and target don't point to the same object, copy all elements of mat1 to target
	if(mat1 != target)
		for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
			*(target->coef + i) = *(mat1->coef + i);
			*(target->expon + i) = *(mat1->expon + i);
		}

	// to multiply numbers in scientific notation, we multiply the coefficients and add the exponents.
	*(target->expon + rowind + colind*mat1->nrow) = *(mat1->expon + rowind + colind*mat1->nrow) + *(mat2->expon);
	*(target->coef + rowind + colind*mat1->nrow) = *(mat1->coef + rowind + colind*mat1->nrow) * *(mat2->coef); 

	// although the sum of the exponents is the approximately correct exponent,
	// the multiplication of coefficients may have thrown things off
	adjust_expon_scinot(target->nrow * target->ncol, target->coef, target->expon);
}


void mult_one_to_one_scinot(scinot_matrix *mat1, int rowind1, int colind1, scinot_matrix *mat2, int rowind2, int colind2, scinot_matrix *target) {
/* multiply element (rowind1, colind1) of mat1 by element (rowind2, colind2) of mat2
 * it is assumed that target has the same dimensions as mat1,
 * rowind1 < mat1->nrow and colind1 < mat1->ncol, and rowind2 < mat2->nrow and colind2 < mat2->ncol
 */
	int i;

	// if mat1 and target don't point to the same object, copy all elements of mat1 to target
	if(mat1 != target)
		for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
			*(target->coef + i) = *(mat1->coef + i);
			*(target->expon + i) = *(mat1->expon + i);
		}

	// to multiply numbers in scientific notation, we multiply the coefficients and add the exponents.
	*(target->expon + rowind1 + colind1*mat1->nrow) = *(mat1->expon + rowind1 + colind1*mat1->nrow) + *(mat2->expon + rowind2 + colind2*mat2->nrow);
	*(target->coef + rowind1 + colind1*mat1->nrow) = *(mat1->coef + rowind1 + colind1*mat1->nrow) * *(mat2->coef + rowind2 + colind2*mat2->nrow); 

	// although the sum of the exponents is the approximately correct exponent,
	// the multiplication of coefficients may have thrown things off
	adjust_expon_scinot(target->nrow * target->ncol, target->coef, target->expon);
}


void mult_to_all_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target) {
/* multiply each element of mat1 by mat2
 * it is assumed that mat1 has the same dimension as target and mat2 is 1 by 1.
 */
	int i;

	// to multiply numbers in scientific notation, we multiply the coefficients and add the exponents.
	for(i = 0; i < mat1->nrow * mat1->ncol; i++) {
		*(target->expon + i) = *(mat1->expon + i) + *(mat2->expon);
		*(target->coef + i) = *(mat1->coef + i) * *(mat2->coef);
	}

	// although the sum of the exponents is the approximately correct exponent,
	// the multiplication of coefficients may have thrown things off
	adjust_expon_scinot(target->nrow * target->ncol, target->coef, target->expon);
}


void matrix_mult_scinot(scinot_matrix *mat1, scinot_matrix *mat2, scinot_matrix *target) {
// calculate the matrix product mat1 mat2 in scientific notation, result stored in target
// for now, it is assumed that all dimensions are correct
// target may be the same as mat1 or mat2 --
// we store the result of the calculation in a temporary variable until the end

	int i, j, k;
	scinot_matrix tempscinot1 = alloc_scinot(1, 1), tempscinotvec1 = alloc_scinot(1, mat1->ncol), tempscinotvec2 = alloc_scinot(1, mat1->ncol), tempres = alloc_scinot(mat1->nrow, mat2->ncol);
//	Rprintf("multiplying matrices of dimensions %i by %i and %i by %i, storing in %i by %i\n", mat1->nrow, mat1->ncol, mat2->nrow, mat2->ncol, target->nrow, target->ncol);

	for(i = 0; i < mat1->nrow; i++) {
		// tempscinotvec1 = row i of mat1
		for(k = 0; k < mat1->ncol; k++) {
			*(tempscinotvec1.coef + k) = *(mat1->coef + i + mat1->nrow * k);
			*(tempscinotvec1.expon + k) = *(mat1->expon + i + mat1->nrow * k);
		}

		for(j = 0; j < mat2->ncol; j++) {
			// tempscinotvec2 = column j of mat2
			for(k = 0; k < mat1->ncol; k++) {
				*(tempscinotvec2.coef + k) = *(mat2->coef + k + mat2->nrow * j);
				*(tempscinotvec2.expon + k) = *(mat2->expon + k + mat2->nrow * j);
			}

			// multiply row i of mat1 by column j of mat2 and sum
			mult_scinot(&tempscinotvec1, &tempscinotvec2, &tempscinotvec2);
			sum_all_entries_scinot(&tempscinotvec2, &tempscinot1);

			// insert to the appropriate location in the result
			insert_scinot(&tempres, &tempscinot1, i, j);
		}
	}

	copy_scinot(&tempres, target);
	
	free_scinot(&tempscinot1);
	free_scinot(&tempscinotvec1);
	free_scinot(&tempscinotvec2);
	free_scinot(&tempres);
}

void log_sum_all_entries_exp_scinot(scinot_matrix *source, scinot_matrix *target) {
/* compute log(exp(a_1) + ... + exp(a_N)), where a_1, ..., a_N are the N elements of the source matrix
 * do this by computing log( sum(a_i)/max(a_i) ) + log(max(a_i))
 */
/*	int i, largestind = 0, len = source->nrow * source->ncol;
	double tempdbl1, tempdbl2;
	scinot_matrix tempscinot1 = alloc_scinot(1, 1), tempscinotDup = alloc_scinot(source->nrow, source->ncol);

	// find the index corresponding to the maximum value
	// easiest way I could think of to compare them was to compute their logarithms
	tempdbl1 = sign(*(source->coef)) * (log10(fabs(*(source->coef))) + *(source->expon));
	for(i = 1; i < len; i++) {
			tempdbl2 = sign(*(source->coef + i)) * (log10(fabs(*(source->coef + i))) + *(source->expon + i));
		if(tempdbl2 > tempdbl1) {
			largestind = i;
			tempdbl1 = tempdbl2;
		}
	}

	// set target = -1 * max(a_i)
	*(target->coef) = -1 * *(source->coef + largestind);
	*(target->expon) = *(source->expon + largestind);

	// subtract max(a_i) from all entries of source
	add_to_all_scinot(source, target, &tempscinotDup);

	// exponentiate, add all entries, and compute the logarithm
	exp_scinot(&tempscinotDup, &tempscinotDup);
	sum_all_entries_scinot(&tempscinotDup, &tempscinot1);
	log_scinot(&tempscinot1, &tempscinot1);
	
	// set target = max(a_i) and add to the other term
	*(target->coef) *= -1;
	add_scinot(target, &tempscinot1, target);

	free_scinot(&tempscinot1);
	free_scinot(&tempscinotDup);
*/
	int i, largestind = 0, len = source->nrow * source->ncol;
	double tempdbl1, tempdbl2;
	
	scinot_matrix tempscinot1 = alloc_scinot(1, 1), tempscinotDup = alloc_scinot(source->nrow, source->ncol);

	// find the index corresponding to the maximum value
	// easiest way I could think of to compare them was to compute their logarithms
	tempdbl1 = sign(*(source->coef)) * (log10(fabs(*(source->coef))) + *(source->expon));
	for(i = 1; i < len; i++) {
			tempdbl2 = sign(*(source->coef + i)) * (log10(fabs(*(source->coef + i))) + *(source->expon + i));
		if(tempdbl2 > tempdbl1) {
			largestind = i;
			tempdbl1 = tempdbl2;
		}
	}

	// set target = -1 * max(a_i)
	*(target->coef) = -1 * *(source->coef + largestind);
	*(target->expon) = *(source->expon + largestind);

	// subtract max(a_i) from all entries of source
	add_to_all_scinot(source, target, &tempscinotDup);

	// exponentiate, add all entries, and compute the logarithm
	exp_scinot(&tempscinotDup, &tempscinotDup);
	sum_all_entries_scinot(&tempscinotDup, &tempscinot1);
	log_scinot(&tempscinot1, &tempscinot1);
	
	// set target = max(a_i) and add to the other term
	*(target->coef) *= -1;
	add_scinot(target, &tempscinot1, target);

	free_scinot(&tempscinot1);
	free_scinot(&tempscinotDup);
}


void exp_scinot(scinot_matrix *source, scinot_matrix *target) {
/* calculates exp(source) in sci notation and stores the result in target
 * assumes that source and target are the same dimensions
 */
	int i, len = source->nrow * source->ncol;
//	double tempvec[len];
	double temp = 0.0;
//	double log10 = log(10);
//	double log_10e = log10(exp(1));
	double log_10e = M_LOG10E;

	// convert each element of logx individually
	for(i = 0; i < len; i++) {
		// temp = log10(x_i)
//		Rprintf("\nlog(x): %lf", *(logx + i));
		if(*(source->expon + i) <= -14 || (*(source->expon + i) == 0 && fabs(*(source->coef + i)) < DBL_EPSILON)) {
//			Rprintf("coef = %lf, DBL_EPSILON = %e", *(source->coef + i), DBL_EPSILON);
			*(target->coef + i) = 0.1;
			*(target->expon + i) = 1;
		} else {
			temp = (log_10e * *(source->coef+i)) * pow(10, *(source->expon + i));
//		Rprintf("\nlog_10(x): %lf", temp);
		// the exponent for x_i
			*(target->expon + i) = ceil(temp);
//		Rprintf("\nexpon: %i", *(expon + i));
		// the coefficient for x_i
			*(target->coef + i) = pow(10, temp - *(target->expon + i));
		}
	}
}

void log_scinot(scinot_matrix *source, scinot_matrix *target) {
/* calculates log(source) in sci notation and stores the result in target
 * assumes that source and target are the same dimensions
 */
	int i, len=source->nrow * source->ncol;
	double log10 = log(10), signcoef;

	// compute log(source) in a sortof numerically stable way -- could be more careful about the exponent
	for(i = 0; i < len; i++) {
		signcoef = sign(*(source->coef + i));
		*(source->coef + i) *= signcoef;
		*(target->coef + i) = signcoef * (log(*(source->coef + i)) + log10 * *(source->expon + i));
		*(target->expon + i) = 0;
	}

	// adjust the coefficients if necessary
	adjust_expon_scinot(len, target->coef, target->expon);
}


void sqrt_scinot(scinot_matrix *source, scinot_matrix *target) {
/* calculates sqrt(source) in sci notation and stores the result in target
 * assumes that source and target are the same dimensions
 */
	int i, len=source->nrow * source->ncol;

	// set the target exponent equal to the source exponent divided by 2
	// take the square root of each coefficient
	// If the source exponent is not divisible by 2, adjust the representation so that it is.
	for(i = 0; i < len; i++) {
		if(abs(*(source->expon + i)) % 2 == 1) {
			*(target->expon + i) = (*(source->expon + i) - 1)/2;
			*(target->coef + i) = sqrt(*(source->coef + i) * 10);
		} else {
			*(target->expon + i) = *(source->expon + i)/2;
			*(target->coef + i) = sqrt(*(source->coef + i));
		}
	}

	// adjust the coefficients if necessary
	adjust_expon_scinot(len, target->coef, target->expon);
}
