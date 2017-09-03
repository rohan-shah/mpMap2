#ifndef COMPUTE_ALL_EPISTATIC_CHI_SQUARED_PROB_HEADER_GUARD
#define COMPUTE_ALL_EPISTATIC_CHI_SQUARED_PROB_HEADER_GUARD
#include <Rcpp.h>
SEXP computeAllEpistaticChiSquared(SEXP probabilities, SEXP nFounders, SEXP infiniteSelfing, SEXP showProgress);
#endif
