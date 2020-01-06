#ifndef GET_FUNNEL_HEADER_GUARD
#define GET_FUNNEL_HEADER_GUARD
#include <Rcpp.h>
void getFunnel(long line, Rcpp::IntegerVector& mother, Rcpp::IntegerVector& father, int* funnel, int nFounders);
#endif

