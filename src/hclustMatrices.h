#ifndef HCLUST_MATRICES_HEADER_GUARD
#define HCLUST_MATRICES_HEADER_GUARD
#include <Rcpp.h>
SEXP hclustCombinedMatrix(SEXP preClusterResults, SEXP mpcrossRF);
SEXP hclustThetaMatrix(SEXP preClusterResults, SEXP mpcrossRF);
SEXP hclustLodMatrix(SEXP preClusterResults, SEXP mpcrossRF);
#endif
