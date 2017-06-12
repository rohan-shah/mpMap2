#ifndef GENERATE_KEYS_HEADER_GUARD
#define GENERATE_KEYS_HEADER_GUARD
#include <Rcpp.h>
void generateKeys(Rcpp::IntegerMatrix& key, Rcpp::IntegerMatrix& outputKey, int nFounders, bool infiniteSelfing);
#endif

