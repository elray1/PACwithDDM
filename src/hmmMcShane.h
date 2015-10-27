/* 
 ============================================================================
 Name        : hmmMcShane.h
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

SEXP calc_log_McShane_class_probs_given_log_static_class_probs_R_interface(SEXP HMM_Xptr, SEXP log_marginal_class_probsp, SEXP Np, SEXP Tp, SEXP log_static_class_probsp);
