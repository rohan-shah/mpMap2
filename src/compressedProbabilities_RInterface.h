#ifndef COMPRESSED_PROBABILITIES_HEADER_GUARD_R_INTERFACE
#define COMPRESSED_PROBABILITIES_HEADER_GUARD_R_INTERFACE
#include <Rcpp.h>
SEXP compressedProbabilities_RInterface(SEXP nFounders, SEXP r, SEXP nFunnels, SEXP intercrossingGenerations, SEXP selfingGenerations, SEXP infiniteSelfing);
#endif
