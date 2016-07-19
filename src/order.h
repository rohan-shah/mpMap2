#ifndef ORDER_HEADER_GUARD
#define ORDER_HEADER_GUARD
#include <Rcpp.h>
SEXP order(SEXP mpcrossLG, SEXP groupsToOrder, SEXP cool_, SEXP temperatureMin_, SEXP nReps_, SEXP maxMove, SEXP effortMultiplier, SEXP randomStart, SEXP verbose_);
#endif
