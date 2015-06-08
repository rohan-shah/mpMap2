#ifndef RECODE_HETS_AS_NA_HEADER_GUARD
#define RECODE_HETS_AS_NA_HEADER_GUARD
#include <Rcpp.h>
bool replaceHetsWithNA(Rcpp::IntegerMatrix recodedFounders, Rcpp::IntegerMatrix recodedFinals, Rcpp::List recodedHetData);
#endif
