/*
 ============================================================================
 Name        : probDistributions.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : probability distributions
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

#include "probDistributions.h"
#include "utility.h"

void rdirich(double *alpha, int len, double *target) {
	/* INPUT:
	 * alpha is a pointer to the first element of an array of length len containing the parameter vector for the dirichlet distribution to sample from.
	 * len is an integer, the dimension of the dirichlet distribution to sample from.
	 *
	 * OUTPUT:
	 * the vector target is modified to contain the sampled value from the dirichlet distribution.
	 * no value is returned.
	 *
	 * The code is based on the comment on page 477 of Cappe, Robert, and Ryden.
	 */
	int i;
	double temp = 0.0;

	for(i = 0; i < len; i++) {
		*(target + i) = rgamma(*(alpha + i), 1);
	}

	for(i = 0; i < len; i++) {
		temp += *(target + i);
	}

	for(i = 0; i < len; i++) {
		*(target + i) /= temp;
	}
}


double ddirich(double *x, double *alpha, int len) {
	/* INPUT:
	 * x is a pointer to the first element of an array of length len containing an observation from the dirichlet distribution.
	 * alpha is a pointer to the first element of an array of length len containing the parameter vector for the dirichlet distribution.
	 * len is an integer, the dimension of the dirichlet distribution to sample from.
	 *
	 * OUTPUT:
	 * returns the density function of the dirichlet(alpha) distribution evaluated at x.
	 */
	int i;
	double retval = 1.0, sum_alpha = 0.0, beta_alpha;

	for(i = 0; i < len; i++) {
		sum_alpha += *(alpha + i);
	}

	beta_alpha = 1/gammafn(sum_alpha);

	for(i = 0; i < len; i++) {
		beta_alpha *= gammafn(*(alpha + i));
		retval *= pow(*(x + i), *(alpha + i) - 1);
	}

	retval /= beta_alpha;

	return retval;
}


void rMVN(int length, double *mu, double *Sigma, double *target) {
	// Compute the square root of Sigma: Same eigenvectors, sqrt of eigenvalues
	int i, j;

	double tempmat[length*length], tempmat2[length*length];
	double sqrtSigma[length*length];

	int eigres, worklen = 4*length;
	double work[worklen], evals[length];
    char jobz[1], uplo[1];
	jobz[0] = 'V';
	uplo[0] = 'L';

	for(i = 0; i < length*length; i++) {
		tempmat[i] = *(Sigma + i);
	}

	// dsyev computes the eigenvalues and eigenvectors of a symmetric matrix
	// eigenvalues are in ascending order, eigenvectors are placed in tempmat
	F77_CALL(dsyev)(jobz, uplo, &length, (double*)&tempmat, &length, (double*)&evals, (double*)&work, &worklen, &eigres);
	
	for(i = 0; i < length; i++) {
		for(j = 0; j < length; j++) {
			if(i == j) {
				sqrtSigma[i + j*length] = sqrt(evals[i]);
			} else {
				sqrtSigma[i + j*length] = 0;
			}
		}
	}


	// compute the square root of Sigma
	char *transno = "N", *transyes = "T";
	double dblone = 1.0, dblzero = 0.0;
	F77_CALL(dgemm)(transno, transno, &length, &length, &length, &dblone, (double*)&tempmat, &length, (double*)&sqrtSigma, &length, &dblzero, (double*)&tempmat2, &length);
	F77_CALL(dgemm)(transno, transyes, &length, &length, &length, &dblone, (double*)&tempmat2, &length, (double*)&tempmat, &length, &dblzero, (double*)&sqrtSigma, &length);

	// Generate length independent normal(0, 1) random variables
	for(i = 0; i < length; i++) {
		tempmat[i] = rnorm(0, 1);
	}

	// Postmultiply the random variables by sqrt(Sigma).  The result is a random vector with mean 0 and the correct covariance matrix
	int intone = 1;
	F77_CALL(dgemm)(transno, transno, &intone, &length, &length, &dblone, (double*)&tempmat, &intone, (double*)&sqrtSigma, &length, &dblzero, target, &intone);

	// Add the mean vector.  The result is a random vector with the correct mean and covariance
	for(i = 0; i < length; i++) {
		*(target + i) += *(mu + i);
	}
}

double dMVN(double *x, int length, int x_pointer_offset, double *mu, double *Sigma, double *log_det_Sigma, int retlog) {
	/* Based on the implementation in the mvtnorm package in R.
	 *
	 * INPUTS:
	 * x is a pointer to an array with the observed value of a multivariate random variable
	 * length is the dimension of x
	 * x_pointer_offset is the number of slots in memory between consecutive entries in x.
	 *    For instance, if x is a row of a matrix that is stored in memory as a vector in column order, x_pointer_offset should be the number of rows in the matrix.
	 * mu is a pointer to an array with the mean vector of x
	 * Sigma is a pointer to an array with the covariance matrix of x
	 * log_det_Sigma is the logarithm of the determinant of Sigma
	 * log is an indicator of whether the log-density or original scale should be returned
	 *
	 * OUTPUTS:
	 * Returns the value of the pdf of a multivariate normal distribution with mean mu and covariance Sigma at the observed x.
	 */

	int i, intone = 1, lwork = -1, info;
	int ipiv[length];
	double logretval, MD = 0.0, opt_lwork;
	double dmmu[length], temp[length], Sigma_copy[length*length];

	char *uplo = "L";
	
//	Rprintf("In dMVN\n");

	/* A call to dsysv with an lwork argument of -1 sets the first element of the work array to the optimal length of the work array */
	F77_CALL(dsysv)(uplo, &length, &intone, (double*)&Sigma_copy, &length, (int*)&ipiv, (double*)&temp, &length, &opt_lwork, &lwork, &info);
	lwork = (int)opt_lwork;
	double work[lwork];

//	Rprintf("Defined workspace of size %i\n", lwork);

	/* make a copy of the Sigma matrix; dsysv destroys the argument passed to it. */
	for(i = 0; i < length*length; i++) {
		Sigma_copy[i] = *(Sigma + i);
	}

//	Rprintf("Copied Sigma\n");
	
	/* make 2 copies of the vector (X - mu) */
	for(i=0; i < length; i++) {
		temp[i] = dmmu[i] = *(x + i*x_pointer_offset) - *(mu + i);
	}
	
//	Rprintf("Copied X - mu, twice\n");

	/* after this call, temp contains Sigma^-1 (X - mu) */
	F77_CALL(dsysv)(uplo, &length, &intone, (double*)&Sigma_copy, &length, (int*)&ipiv, (double*)&temp, &length, (double*)&work, &lwork, &info);

//	Rprintf("Computation step one\n");
	
	/* compute (X - mu)' Sigma^-1 (X - mu) */
	for(i=0; i < length; i++) {
		MD += dmmu[i] * temp[i];
	}

//	Rprintf("Computation step two\n");
	
	logretval = -0.5*( length * log(2*M_PI) + *log_det_Sigma + MD );
	
//	Rprintf("Got loglik\n");

	if(retlog) {
		return logretval;
	} else {
//		Rprintf("Exponentiating result in C dMVN\n");
		return exp(logretval);
	}
}

SEXP dMVN_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog) {
	SEXP retval = PROTECT(allocVector(REALSXP, *INTEGER(n)));
	
	dMVN_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
	
	UNPROTECT(1);
	return(retval);
}

SEXP dMVN_Prec_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog) {
	SEXP retval = PROTECT(allocVector(REALSXP, *INTEGER(n)));
	
	dMVN_Prec_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
	
	UNPROTECT(1);
	return(retval);
}

SEXP dMVN_DiagPrec_multiobs_R_interface(SEXP x, SEXP n, SEXP d, SEXP mu, SEXP Sigma, SEXP log_det_Sigma, SEXP retlog) {
	SEXP retval = PROTECT(allocVector(REALSXP, *INTEGER(n)));
	
	dMVN_DiagPrec_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
	
	UNPROTECT(1);
	return(retval);
}


void dMVN_multiobs(double *x, int n, int d, double *mu, double *Sigma, double *log_det_Sigma, double *result, int retlog) {
	/* Based on the implementation in the mvtnorm package in R.
	 *
	 * INPUTS:
	 * x is a pointer to an n by d array with the observed value of a multivariate random variable
	 * length is the dimension of x
	 * x_pointer_offset is the number of slots in memory between consecutive entries in x.
	 *    For instance, if x is a row of a matrix that is stored in memory as a vector in column order, x_pointer_offset should be the number of rows in the matrix.
	 * mu is a pointer to an array with the mean vector of x
	 * Sigma is a pointer to an array with the covariance matrix of x
	 * log_det_Sigma is the logarithm of the determinant of Sigma
	 * log is an indicator of whether the log-density or original scale should be returned
	 *
	 * OUTPUTS:
	 * Returns the value of the pdf of a multivariate normal distribution with mean mu and covariance Sigma at the observed x.
	 */

	int i, j, lwork = -1, info, ni = n, di = d;
	int ipiv[d];
	double opt_lwork;
	double Sigma_copy[d * d];
//	Rprintf("n = %i, d = %i\n", n, d);
	double *dmmu = Calloc(n * d, double);
	double *temp = Calloc(n * d, double);

	char *uplo = "L";
	
//	Rprintf("In dMVN\n");

	/* A call to dsysv with an lwork argument of -1 sets the first element of the work array to the optimal length of the work array */
	F77_CALL(dsysv)(uplo, &di, &ni, (double*)&Sigma_copy, &di, (int*)&ipiv, temp, &di, &opt_lwork, &lwork, &info);
	lwork = (int)opt_lwork;
	double work[lwork];

//	Rprintf("Defined workspace of size %i\n", lwork);

	/* make a copy of the Sigma matrix; dsysv destroys the argument passed to it. */
	for(i = 0; i < d * d; i++) {
		Sigma_copy[i] = *(Sigma + i);
	}
	
//	Rprintf("Copied Sigma\n");
	
	/* make 2 copies of the vector (X - mu).  We transpose the order of the elements to get a d by n matrix */
	for(i = 0; i < n; i++) {
		for(j = 0; j < d; j++) {
			*(temp + j + i*d) = *(dmmu + j + i*d) = *(x + i + n*j) - *(mu + j);
//			dmmu[j + i*d] = *(x + i + n*j) - *(mu + j);
		}
	}
	
//	Rprintf("Copied X - mu, twice\n");
	
	/* after this call, dmmu contains Sigma^-1 (X - mu)' */
	F77_CALL(dsysv)(uplo, &d, &n, (double*)&Sigma_copy, &d, (int*)&ipiv, dmmu, &d, (double*)&work, &lwork, &info);
	
//	Rprintf("Computation step one\n");
	
	/* compute (X - mu)' Sigma^-1 (X - mu) */
	double norm_const = d * log(2 * M_PI) + *log_det_Sigma; 
	for(i = 0; i < n; i++) {
		*(result + i) = norm_const;
		for(j = 0; j < d; j++) {
			Rprintf("%lf ", *(dmmu + j + i * d));
			*(result + i) += *(dmmu + j + i * d) * *(temp + j + i * d);
		}
		*(result + i) *= -0.5;
		Rprintf(" -> %lf\n", *(result + i));
	}

//	matrix_multiply((double *)&temp, (double *)&dmmu, result, n, d, n);
	
//	Rprintf("Computation step two\n");
	
//	for(i = 0; i < n; i++) {
//		*(result + i) = -0.5*( d * log(2*M_PI) + *log_det_Sigma + *(result + i) );
//	}
	
//	Rprintf("Got loglik\n");
	
	if(!retlog) {
//		Rprintf("Exponentiating result in C dMVN\n");
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(dmmu);
	Free(temp);
}


void dMVN_Prec_multiobs(double *x, int n, int d, double *mu, double *Prec, double *log_det_Sigma, double *result, int retlog) {
	/* Based on the implementation in the mvtnorm package in R.
	 *
	 * INPUTS:
	 * x is a pointer to an n by d array with the observed value of a multivariate random variable
	 * length is the dimension of x
	 * x_pointer_offset is the number of slots in memory between consecutive entries in x.
	 *    For instance, if x is a row of a matrix that is stored in memory as a vector in column order, x_pointer_offset should be the number of rows in the matrix.
	 * mu is a pointer to an array with the mean vector of x
	 * Sigma is a pointer to an array with the covariance matrix of x
	 * log_det_Sigma is the logarithm of the determinant of Sigma
	 * log is an indicator of whether the log-density or original scale should be returned
	 *
	 * OUTPUTS:
	 * Returns the value of the pdf of a multivariate normal distribution with mean mu and covariance Sigma at the observed x.
	 */

	int i, j;
	double *dmmu = Calloc(n * d, double), *temp = Calloc(n * d, double);
	
//	Rprintf("In dMVN\n");

	/* A call to dsysv with an lwork argument of -1 sets the first element of the work array to the optimal length of the work array */
//	F77_CALL(dsysv)(uplo, &di, &ni, (double*)&Sigma_copy, &di, (int*)&ipiv, temp, &di, &opt_lwork, &lwork, &info);
//	lwork = (int)opt_lwork;
//	double work[lwork];

//	Rprintf("Defined workspace of size %i\n", lwork);

//	/* make a copy of the Sigma matrix; dsysv destroys the argument passed to it. */
//	for(i = 0; i < d * d; i++) {
//		Sigma_copy[i] = *(Sigma + i);
//	}
	
//	Rprintf("Copied Sigma\n");
	
	/* make 2 copies of the vector (X - mu).  We transpose the order of the elements to get a d by n matrix */
	for(i = 0; i < n; i++) {
		for(j = 0; j < d; j++) {
			*(temp + j + i*d) = *(dmmu + j + i*d) = *(x + i + n*j) - *(mu + j);
//			dmmu[j + i*d] = *(x + i + n*j) - *(mu + j);
		}
	}
	
//	Rprintf("Copied X - mu, twice\n");
	
	/* after this call, dmmu contains Sigma^-1 (X - mu)' */
	matrix_multiply(Prec, dmmu, dmmu, d, d, n);
	//	F77_CALL(dsysv)(uplo, &d, &n, (double*)&Sigma_copy, &d, (int*)&ipiv, dmmu, &d, (double*)&work, &lwork, &info);
	
//	Rprintf("Computation step one\n");
	
	/* compute (X - mu)' Sigma^-1 (X - mu) */
	double norm_const = d * log(2 * M_PI) + *log_det_Sigma; 
	for(i = 0; i < n; i++) {
		*(result + i) = norm_const;
		for(j = 0; j < d; j++) {
			*(result + i) += *(dmmu + j + i * d) * *(temp + j + i * d);
		}
		*(result + i) *= -0.5;
	}
	
//	Rprintf("Computation step two\n");
	
//	for(i = 0; i < n; i++) {
//		*(result + i) = -0.5*( d * log(2*M_PI) + *log_det_Sigma + *(result + i) );
//	}
	
//	Rprintf("Got loglik\n");
	
	if(!retlog) {
//		Rprintf("Exponentiating result in C dMVN\n");
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(dmmu);
	Free(temp);
}


void dMVN_DiagPrec_multiobs(double *x, int n, int d, double *mu, double *Prec, double *log_det_Sigma, double *result, int retlog) {
	/* Based on the implementation in the mvtnorm package in R.
	 *
	 * INPUTS:
	 * x is a pointer to an n by d array with the observed value of a multivariate random variable
	 * length is the dimension of x
	 * x_pointer_offset is the number of slots in memory between consecutive entries in x.
	 *    For instance, if x is a row of a matrix that is stored in memory as a vector in column order, x_pointer_offset should be the number of rows in the matrix.
	 * mu is a pointer to an array with the mean vector of x
	 * Sigma is a pointer to an array with the covariance matrix of x
	 * log_det_Sigma is the logarithm of the determinant of Sigma
	 * log is an indicator of whether the log-density or original scale should be returned
	 *
	 * OUTPUTS:
	 * Returns the value of the pdf of a multivariate normal distribution with mean mu and covariance Sigma at the observed x.
	 */

	int i, j;
	double *dmmu = Calloc(n * d, double), *temp = Calloc(n * d, double);
	
//	Rprintf("In dMVN\n");

	/* A call to dsysv with an lwork argument of -1 sets the first element of the work array to the optimal length of the work array */
//	F77_CALL(dsysv)(uplo, &di, &ni, (double*)&Sigma_copy, &di, (int*)&ipiv, temp, &di, &opt_lwork, &lwork, &info);
//	lwork = (int)opt_lwork;
//	double work[lwork];

//	Rprintf("Defined workspace of size %i\n", lwork);

//	/* make a copy of the Sigma matrix; dsysv destroys the argument passed to it. */
//	for(i = 0; i < d * d; i++) {
//		Sigma_copy[i] = *(Sigma + i);
//	}
	
//	Rprintf("Copied Sigma\n");
	
	/* make 2 copies of the vector (X - mu).  We transpose the order of the elements to get a d by n matrix */
	for(i = 0; i < n; i++) {
		for(j = 0; j < d; j++) {
			*(dmmu + j + i*d) = *(x + i + n*j) - *(mu + j);
//			dmmu[j + i*d] = *(x + i + n*j) - *(mu + j);
		}
	}
	
//	Rprintf("Copied X - mu, twice\n");
	
	/* after this call, dmmu contains Sigma^-1 (X - mu)' */
//	matrix_multiply(Prec, dmmu, dmmu, d, d, n);
	//	F77_CALL(dsysv)(uplo, &d, &n, (double*)&Sigma_copy, &d, (int*)&ipiv, dmmu, &d, (double*)&work, &lwork, &info);
	
//	Rprintf("Computation step one\n");
	
	/* compute (X - mu)' Sigma^-1 (X - mu) */
	double norm_const = d * log(2 * M_PI) + *log_det_Sigma; 
	for(i = 0; i < n; i++) {
		*(result + i) = norm_const;
		for(j = 0; j < d; j++) {
			*(result + i) += *(dmmu + j + i * d) * *(dmmu + j + i * d) * *(Prec + j);
		}
		*(result + i) *= -0.5;
	}
	
//	Rprintf("Computation step two\n");
	
//	for(i = 0; i < n; i++) {
//		*(result + i) = -0.5*( d * log(2*M_PI) + *log_det_Sigma + *(result + i) );
//	}
	
//	Rprintf("Got loglik\n");
	
	if(!retlog) {
//		Rprintf("Exponentiating result in C dMVN\n");
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(dmmu);
	Free(temp);
}


SEXP dGMM_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Sigmasp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog) {
//	Rprintf("made it into dGMM\n");
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp);
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Sigma, *log_det_Sigma;
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
//	Rprintf("M = %i\n", M);
	
	if(M > 0) {
		Rprintf("model 0 of %i\n", M);
		log_rho = log(*(rhos));
		Rprintf("log rho = %lf\n", log_rho);
		mu = REAL(VECTOR_ELT(musp, 0));
		Sigma = REAL(VECTOR_ELT(Sigmasp, 0));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, 0));
		
		dMVN_multiobs(x, n, d, mu, Sigma, log_det_Sigma, result, 1);
		for(i = 0; i < n; i++) {
			Rprintf("dmvn obs i = %i: %lf, result = ", i, *(result + i));
			*(result + i) += log_rho;
			Rprintf("%lf\n", *(result + i));
		}
	}
	
	for(m = 1; m < M; m++) {
		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		Rprintf("log rho = %lf\n", log_rho);
		mu = REAL(VECTOR_ELT(musp, m));
		Sigma = REAL(VECTOR_ELT(Sigmasp, m));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, m));
		
		dMVN_multiobs(x, n, d, mu, Sigma, log_det_Sigma, tempdbl, 1);
		
		for(i = 0; i < n; i++) {
			Rprintf("dmvn obs i = %i: %lf, combined result = ", i, *(tempdbl + i));
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
			Rprintf("%lf\n", *(result + i));
		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	
	return(resultsexp);
}

SEXP dGMM_sameSigma_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Sigmap, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog) {
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp), intone = 1;
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Sigma = REAL(Sigmap), *log_det_Sigma = REAL(log_det_Sigmap);
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
	
	if(M > 0) {
//		Rprintf("model %i of %i\n", 0, M);
		log_rho = log(*(rhos));
		mu = REAL(VECTOR_ELT(musp, 0));
		
//		dMVN_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
		dMVN_multiobs(x, n, d, mu, Sigma, log_det_Sigma, result, intone);
		for(i = 0; i < n; i++) {
			*(result + i) += log_rho;
		}
//		for(i = 0; i < n; i++) {
//			Rprintf("i = %i\n", i);
//			*(result + i) = dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho;
//		}
	}
	
	for(m = 1; m < M; m++) {
//		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		mu = REAL(VECTOR_ELT(musp, m));
		
		dMVN_multiobs(x, n, d, mu, Sigma, log_det_Sigma, tempdbl, 1);
		
		for(i = 0; i < n; i++) {
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
		}
//		for(i = 0; i < n; i++) {
//			*(result + i) = logspace_add(*(result + i),
//					dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho);
//		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	return(resultsexp);
}


SEXP dGMM_Prec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precsp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog) {
//	Rprintf("made it into dGMM\n");
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp);
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Prec, *log_det_Sigma;
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
//	Rprintf("M = %i\n", M);
	
	if(M > 0) {
//		Rprintf("model 0 of %i\n", M);
		log_rho = log(*(rhos));
		mu = REAL(VECTOR_ELT(musp, 0));
		Prec = REAL(VECTOR_ELT(Precsp, 0));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, 0));
		
		dMVN_Prec_multiobs(x, n, d, mu, Prec, log_det_Sigma, result, 1);
		for(i = 0; i < n; i++) {
			*(result + i) += log_rho;
		}
	}
	
	for(m = 1; m < M; m++) {
//		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		mu = REAL(VECTOR_ELT(musp, m));
		Prec = REAL(VECTOR_ELT(Precsp, m));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, m));
		
		dMVN_Prec_multiobs(x, n, d, mu, Prec, log_det_Sigma, tempdbl, 1);
		
		for(i = 0; i < n; i++) {
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	
	return(resultsexp);
}

SEXP dGMM_samePrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precp, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog) {
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp), intone = 1;
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Prec = REAL(Precp), *log_det_Sigma = REAL(log_det_Sigmap);
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
	
	if(M > 0) {
//		Rprintf("model %i of %i\n", 0, M);
		log_rho = log(*(rhos));
		mu = REAL(VECTOR_ELT(musp, 0));
		
//		dMVN_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
		dMVN_Prec_multiobs(x, n, d, mu, Prec, log_det_Sigma, result, intone);
		for(i = 0; i < n; i++) {
			*(result + i) += log_rho;
		}
//		for(i = 0; i < n; i++) {
//			Rprintf("i = %i\n", i);
//			*(result + i) = dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho;
//		}
	}
	
	for(m = 1; m < M; m++) {
//		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		mu = REAL(VECTOR_ELT(musp, m));
		
		dMVN_Prec_multiobs(x, n, d, mu, Prec, log_det_Sigma, tempdbl, 1);
		
		for(i = 0; i < n; i++) {
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
		}
//		for(i = 0; i < n; i++) {
//			*(result + i) = logspace_add(*(result + i),
//					dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho);
//		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	return(resultsexp);
}


SEXP dGMM_DiagPrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precsp, SEXP log_det_Sigmasp, SEXP Mp, SEXP retlog) {
//	Rprintf("made it into dGMM\n");
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp);
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Prec, *log_det_Sigma;
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
//	Rprintf("M = %i\n", M);
	
	if(M > 0) {
//		Rprintf("model 0 of %i\n", M);
		log_rho = log(*(rhos));
		mu = REAL(VECTOR_ELT(musp, 0));
		Prec = REAL(VECTOR_ELT(Precsp, 0));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, 0));
		
		dMVN_DiagPrec_multiobs(x, n, d, mu, Prec, log_det_Sigma, result, 1);
		for(i = 0; i < n; i++) {
			*(result + i) += log_rho;
		}
	}
	
	for(m = 1; m < M; m++) {
//		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		mu = REAL(VECTOR_ELT(musp, m));
		Prec = REAL(VECTOR_ELT(Precsp, m));
		log_det_Sigma = REAL(VECTOR_ELT(log_det_Sigmasp, m));
		
		dMVN_DiagPrec_multiobs(x, n, d, mu, Prec, log_det_Sigma, tempdbl, 1);
		
		for(i = 0; i < n; i++) {
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	
	return(resultsexp);
}

SEXP dGMM_sameDiagPrec_R_Interface(SEXP Xp, SEXP np, SEXP dp, SEXP rhosp, SEXP musp, SEXP Precp, SEXP log_det_Sigmap, SEXP Mp, SEXP retlog) {
	int i, m, n = *INTEGER(np), d = *INTEGER(dp), M = *INTEGER(Mp), intone = 1;
	double log_rho;
	double *x = REAL(Xp);
	double *rhos = REAL(rhosp), *mu, *Prec = REAL(Precp), *log_det_Sigma = REAL(log_det_Sigmap);
	
	SEXP resultsexp = PROTECT(allocVector(REALSXP, n));
	double *result = REAL(resultsexp);
	double *tempdbl = Calloc(n, double);
	
	if(M > 0) {
//		Rprintf("model %i of %i\n", 0, M);
		log_rho = log(*(rhos));
		mu = REAL(VECTOR_ELT(musp, 0));
		
//		dMVN_multiobs(REAL(x), *INTEGER(n), *INTEGER(d), REAL(mu), REAL(Sigma), REAL(log_det_Sigma), REAL(retval), *INTEGER(retlog));
		dMVN_DiagPrec_multiobs(x, n, d, mu, Prec, log_det_Sigma, result, intone);
//		Rprintf("log rho = %lf\n", log_rho);
//		for(i = 0; i < 100; i++)
//			Rprintf("%lf, ", *(tempdbl + i));
//		Rprintf("\n");
		for(i = 0; i < n; i++) {
			*(result + i) += log_rho;
		}
//		for(i = 0; i < n; i++) {
//			Rprintf("i = %i\n", i);
//			*(result + i) = dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho;
//		}
	}
	
	for(m = 1; m < M; m++) {
//		Rprintf("model %i of %i\n", m, M);
		log_rho = log(*(rhos + m));
		mu = REAL(VECTOR_ELT(musp, m));
		
		dMVN_DiagPrec_multiobs(x, n, d, mu, Prec, log_det_Sigma, tempdbl, 1);
//		Rprintf("log rho = %lf\n", log_rho);
//		for(i = 0; i < 100; i++)
//			Rprintf("%lf, ", *(tempdbl + i));
//		Rprintf("\n");
		
		for(i = 0; i < n; i++) {
			*(result + i) = logspace_add(*(result + i), *(tempdbl + i) + log_rho);
		}
//		for(i = 0; i < n; i++) {
//			*(result + i) = logspace_add(*(result + i),
//					dMVN(x + i, d, n, mu, Sigma, log_det_Sigma, 1) + log_rho);
//		}
	}
	
	if(!*INTEGER(retlog)) {
		for(i = 0; i < n; i++) {
			*(result + i) = exp(*(result + i));
		}
	}
	
	Free(tempdbl);
	UNPROTECT(1);
	return(resultsexp);
}


int uneq_prob_sample(double *prop_prob_vec, int length) {
	double tempdbl = 0.0, prob_sums[length];
	int i;

	/* calculate the cumulative probabilities, in the vector prob_sums */
	for(i = 0; i < length; i++) {
		tempdbl += *(prop_prob_vec + i);
	}

	prob_sums[0] = *(prop_prob_vec + 0)/tempdbl;
	for(i = 1; i < length; i++) {
		prob_sums[i] = *(prop_prob_vec + i)/tempdbl + prob_sums[i-1];
	}

	/* draw the sample */
	tempdbl = runif(0, 1);
	for(i = 0; i < length; i++) {
		if(tempdbl < prob_sums[i]) {
			break;
		}
	}

	return i;
}

int runif_int(int min, int max) {
	if(min == max)
		return min;

	return min + (int)( ((double)(max - min + 1)) * runif(0, 1) );
}

void sample_int_without_replace_FisherYates(int *values, int K, int *result, int N) {
	// sample N integers without replacement from a set of K integers
	// values is an array of length K containing the integers to sample from
	// result is an array of length N that will be used to store the results
	int max = K - 1, r;

	while(N >= 1) {
//Rprintf("N = %i\n", N);
		// choose an index in value between 0 and max (inclusive).
		r = runif_int(0, max);

		// store the chosen value in the result vector
		*(result + N - 1) = *(values + r);

		// swap the chosen value with the value in position max
		*(values + r) = *(values + max);
		*(values + max) = *(result + N - 1);

		// decrement N and max
		N--;
		max--;
	}
}

void sample_int_without_replace(int min, int max, int N, int *result) {
	if(N >= 1) {
		int i, j, K, found;
		
		K = max - min + 1;

		int *selected_table = Calloc(N, int);

		if(N == K) {
			for(i = 0; i < N; i++)
				*(result + i) = i;
		} else if(N < K) {
			if(N > 1)
				Rprintf("K = %i, N = %i, result is: ", K, N);
			for(i = 0; i < N; i++) {
				// which of the remaining K empty slots will we fill?
				*(selected_table + i) = runif_int(0, K);

				// if this value was selected previously, map it to the corresponding unchosen value.
				found = 0;
				for(j = i - 1; j >= 0; j--) {
					if(*(selected_table + j) == *(selected_table + i)) {
						*(result + i) = min + K - j;
					}
				}

				// otherwise, store this value
				if(!found)
					*(result + i) = min + *(selected_table + i);

				if(N > 1)
					Rprintf("%i, ", *(result + i));
			}
			if(N > 1)
				Rprintf("\n");
		}
	}
}
