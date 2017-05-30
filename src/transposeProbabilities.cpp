#include "transposeProbabilities.h"
#include <iterator>
#include <vector>
#include <string>
SEXP transposeProbabilities(SEXP probabilities_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData;
	try
	{
		geneticData = geneticData_sexp;
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData must be an S4 object");
	}

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@finals must be an integer matrix");
	}
	int nFinals = finals.nrow();
	
	Rcpp::IntegerMatrix founders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@founders must be an integer matrix");
	}

	Rcpp::S4 probabilities;
	try
	{
		probabilities = Rcpp::as<Rcpp::S4>(geneticData.slot("probabilities"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@probabilities must be an S4 object");
	}

	Rcpp::NumericMatrix data;
	try
	{
		data = Rcpp::as<Rcpp::NumericMatrix>(probabilities.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@probabilities@data must be a numeric matrix");
	}

	Rcpp::List map;
	try
	{
		map = Rcpp::as<Rcpp::List>(probabilities.slot("map"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@probabilities@map must be a list");
	}

	Rcpp::IntegerMatrix key;
	try
	{
		key = Rcpp::as<Rcpp::IntegerMatrix>(probabilities.slot("key"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@probabilities@key must be an integer matrix");
	}
	std::vector<int> keyValues;
	keyValues.insert(keyValues.begin(), key.column(2).begin(), key.column(2).end());
	std::sort(keyValues.begin(), keyValues.end());
	keyValues.erase(std::unique(keyValues.begin(), keyValues.end()), keyValues.end());
	int nAlleles = keyValues.size();

	std::vector<std::string> positionNames;
	for(Rcpp::List::iterator i = map.begin(); i != map.end(); i++)
	{
		Rcpp::CharacterVector names = Rcpp::as<Rcpp::NumericVector>(*i).names();
		std::transform(names.begin(), names.end(), std::back_inserter(positionNames), &Rcpp::as<std::string>);
	}

	if(data.ncol() != positionNames.size() || data.nrow() != nAlleles * nFinals)
	{
		throw std::runtime_error("Input geneticData@probabilities@data had the wrong dimensions");
	}

	Rcpp::NumericMatrix results(nFinals, positionNames.size()*nAlleles);
	Rcpp::CharacterVector newColumnNames(positionNames.size()*nAlleles);
	for(int alleleCounter = 0; alleleCounter < nAlleles; alleleCounter++)
	{
		for(int positionCounter = 0; positionCounter < positionNames.size(); positionCounter++)
		{
			std::stringstream ss;
			ss << positionNames[positionCounter] << " - " << (alleleCounter+1);
			newColumnNames(positionCounter * nAlleles + alleleCounter) = ss.str();
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				results(finalCounter, positionCounter*nAlleles+ alleleCounter) = data(finalCounter*nAlleles+ alleleCounter, positionCounter);
			}
		}
	}

	Rcpp::rownames(results) = Rcpp::rownames(finals);
	Rcpp::colnames(results) = newColumnNames;
	return(results);
END_RCPP
}

