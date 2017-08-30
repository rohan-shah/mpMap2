#ifndef MPMAP2_ARSA_RAW_R_ENTRY_POINT_HEADER_GUARD
#define MPMAP2_ARSA_RAW_R_ENTRY_POINT_HEADER_GUARD
#include "Rcpp.h"
#include <functional>
SEXP arsaRawREntryPoint(SEXP n_, SEXP rawDist_, SEXP levels_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_, SEXP maxMove_sexp, SEXP effortMultiplier_sexp, SEXP randomStart_sexp);
#endif

