#ifndef ALLELE_DATA_ERRORS_HEADER_GUARD
#define ALLELE_DATA_ERRORS_HEADER_GUARD
#include <Rcpp.h>
RcppExport SEXP alleleDataErrors(SEXP Robject);
RcppExport SEXP listCodingErrors(SEXP founders, SEXP finals, SEXP hetData);
#endif