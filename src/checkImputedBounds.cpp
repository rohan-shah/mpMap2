#include "checkImputedBounds.h"
SEXP checkImputedBounds(SEXP matrix_sexp, SEXP upperBound_sexp)
{
BEGIN_RCPP
	Rcpp::IntegerVector matrix;
	try
	{
		matrix = Rcpp::as<Rcpp::IntegerVector>(matrix_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input matrix must be an integer matrix");
	}

	int upperBound;
	try
	{
		upperBound = Rcpp::as<int>(upperBound_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input upperBound must be an integer");
	}
	for(int i = 0; i < matrix.size(); i++)
	{
		if(matrix(i) < 1 || matrix(i) > upperBound)
		{
			return Rcpp::wrap(false);
		}
	}
	return Rcpp::wrap(true);
END_RCPP
}
