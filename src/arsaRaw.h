#ifndef MPMAP2_ARSA_HEADER_GUARD
#define MPMAP2_ARSA_HEADER_GUARD
#include "Rcpp.h"
#include <functional>
SEXP arsaRaw(SEXP n_, SEXP rawDist_, SEXP levels_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_);
void arsaRaw(long n, Rbyte* rawDist, std::vector<double>& levels, double cool, double temperatureMin, long nReps, std::vector<int>& permutation, std::function<void(unsigned long,unsigned long)> progressFunction);
#endif

