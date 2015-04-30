#ifndef RECODE_FOUNDERS_AND_FINALS_AND_HET_DATA_HEADER_GUARD
#define RECODE_FOUNDERS_AND_FINALS_AND_HET_DATA_HEADER_GUARD
#include <Rcpp.h>
void recodeFoundersFinalsHets(Rcpp::IntegerMatrix& recodedFounders, Rcpp::IntegerMatrix& recodedFinals, Rcpp::IntegerMatrix& founders, Rcpp::IntegerMatrix& finals, unsigned int& maxAlleles);
#endif