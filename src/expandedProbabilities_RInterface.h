#ifndef EXPANDED_PROBABILITIES_INFINITE_HEADER_GUARD_R_INTERFACE
#define EXPANDED_PROBABILITIES_INFINITE_HEADER_GUARD_R_INTERFACE
#include <Rcpp.h>
SEXP expandedProbabilitiesInfinite_RInterface(SEXP nFounders, SEXP r, SEXP nFunnels, SEXP intercrossingGenerations);
SEXP expandedProbabilitiesFinite_RInterface(SEXP nFounders_sexp, SEXP r_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp, SEXP selfingGenerations_sexp, SEXP phased_sexp);
#endif
