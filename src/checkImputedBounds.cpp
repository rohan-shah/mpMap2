#include "checkImputedBounds.h"
SEXP checkImputedBounds(SEXP imputed_sexp)
{
BEGIN_RCPP
	Rcpp::S4 imputed;
	try
	{
		imputed = Rcpp::as<Rcpp::S4>(imputed_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input imputed must be an S4 object");
	}

	Rcpp::IntegerVector dataMatrix;
	try
	{
		dataMatrix = Rcpp::as<Rcpp::IntegerVector>(imputed.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Input imputed@data must be an integer matrix");
	}

	Rcpp::IntegerMatrix keyMatrix;
	try
	{
		keyMatrix = Rcpp::as<Rcpp::IntegerMatrix>(imputed.slot("key"));
	}
	catch(...)
	{
		throw std::runtime_error("Input imputed@key must be an integer matrix");
	}

	if(keyMatrix.ncol() != 3)
	{
		throw std::runtime_error("Input imputed@key must have three columns");
	}

	std::vector<int> encodings;
	for(int i = 0; i < keyMatrix.nrow(); i++)
	{
		encodings.push_back(keyMatrix(i, 2));
	}
	std::sort(encodings.begin(), encodings.end());
	encodings.erase(std::unique(encodings.begin(), encodings.end()), encodings.end());
	if(encodings[0] != 1 || *encodings.rbegin() != (int)encodings.size())
	{
		throw std::runtime_error("Encodings in imputed@key must be consecutive integers");
	}

	int upperBound = (int)encodings.size();
	for(int i = 0; i < dataMatrix.size(); i++)
	{
		if(dataMatrix(i) < 1 || dataMatrix(i) > upperBound)
		{
			return Rcpp::wrap(false);
		}
	}
	return Rcpp::wrap(true);
END_RCPP
}
