#include "alleleDataErrors.h"
#include <vector>
#include <string>
RcppExport SEXP alleleDataErrors(SEXP Robject)
{
	BEGIN_RCPP
		Rcpp::S4 mpcross = Robject;
		Rcpp::IntegerMatrix founders = mpcross.slot("founders"), finals = mpcross.slot("finals");
		Rcpp::S4 hetData = mpcross.slot("hetData");

		Rcpp::List codingErrors = listCodingErrors(founders, finals, hetData);

		std::vector<std::string> errors;
		errors.push_back("Typical error");
		return Rcpp::wrap(errors);
	END_RCPP
}
RcppExport SEXP listCodingErrors(SEXP _founders, SEXP _finals, SEXP _hetData)
{
	BEGIN_RCPP
		Rcpp::IntegerMatrix founders = _founders, finals = _finals;
		Rcpp::List hetData = _hetData;

		std::vector<int> _founderErrors, _finalErrors;
		Rcpp::IntegerVector founderErrors = Rcpp::wrap(_founderErrors);
		Rcpp::IntegerVector finalErrors = Rcpp::wrap(_finalErrors);
		return Rcpp::List::create(Rcpp::Named("founders") = founderErrors, Rcpp::Named("finals") = finalErrors);
	END_RCPP
}