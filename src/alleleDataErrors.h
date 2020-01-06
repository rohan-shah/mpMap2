#ifndef ALLELE_DATA_ERRORS_HEADER_GUARD
#define ALLELE_DATA_ERRORS_HEADER_GUARD
#include <Rcpp.h>
RcppExport SEXP alleleDataErrors(SEXP Robject, SEXP Rlimit);
RcppExport void codingErrorsToStrings(Rcpp::List codingErrors, std::vector<std::string>& codingErrorsAsStrings, Rcpp::IntegerMatrix finals, Rcpp::List hetData, int limit);
RcppExport SEXP listCodingErrors(SEXP founders, SEXP finals, SEXP hetData);
#endif

