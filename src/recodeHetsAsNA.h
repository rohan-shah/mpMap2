#ifndef RECODE_HETS_AS_NA_HEADER_GUARD
#define RECODE_HETS_AS_NA_HEADER_GUARD
#include <Rcpp.h>
void replaceHetsWithNA(Rcpp::IntegerMatrix recodedFounders, Rcpp::IntegerMatrix recodedFinals, Rcpp::List recodedHetData, bool& hasHets, bool& hasHetEncodings);
#endif
