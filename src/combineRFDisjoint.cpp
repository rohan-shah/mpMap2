#include "combineRFDisjoint.h"
RcppExport SEXP combineRFDisjoint(SEXP rf1_sexp, SEXP rf2_sexp)
{
BEGIN_RCPP
	Rcpp::S4 rf1;
	try
	{
		rf1 = Rcpp::as<Rcpp::S4>(rf1_sexp);
	}
	catch(...)
	{
	}

	Rcpp::S4 rf2;
	try
	{
		rf2 = Rcpp::as<Rcpp::S4>(rf2_sexp);
	}
	catch(...)
	{
	}

	Rcpp::S4 theta1 = Rcpp::as<Rcpp::S4>(rf1.slot("theta"));
	Rcpp::S4 theta2 = Rcpp::as<Rcpp::S4>(rf2.slot("theta"));

	Rcpp::CharacterVector markers1 = Rcpp::as<Rcpp::CharacterVector>(theta1.slot("markers"));
	Rcpp::CharacterVector markers2 = Rcpp::as<Rcpp::CharacterVector>(theta2.slot("markers"));

	Rcpp::NumericVector levels1 = Rcpp::as<Rcpp::NumericVector>(theta1.slot("levels"));
	Rcpp::NumericVector levels2 = Rcpp::as<Rcpp::NumericVector>(theta2.slot("levels"));

	Rcpp::RawVector raw1 = Rcpp::as<Rcpp::RawVector>(theta1.slot("data"));
	Rcpp::RawVector raw2 = Rcpp::as<Rcpp::RawVector>(theta2.slot("data"));

	if(levels1.size() != levels2.size() || !std::equal(levels1.begin(), levels1.end(), levels2.begin()))
	{
		throw std::runtime_error("Levels must be identical for slot theta in both objects");
	}
	
	std::size_t nMarkers1 = (std::size_t)markers1.size();
	std::size_t nMarkers2 = (std::size_t)markers2.size();
	std::size_t newMarkers = (std::size_t)nMarkers1 + nMarkers2;

	Rcpp::RawVector combinedRaw(newMarkers * (newMarkers+1ULL)/2ULL);
	std::fill(combinedRaw.begin(), combinedRaw.end(), -1);
	for(std::size_t columnCounter = 0; columnCounter < nMarkers1; columnCounter++)
	{
		for(std::size_t rowCounter = 0; rowCounter <= columnCounter; rowCounter++)
		{
			combinedRaw(columnCounter*(columnCounter+1ULL)/2ULL + rowCounter) = raw1(columnCounter*(columnCounter+1ULL)/2ULL + rowCounter); 
		}
	}

	for(std::size_t columnCounter = 0; columnCounter < nMarkers2; columnCounter++)
	{
		for(std::size_t rowCounter = 0; rowCounter <= columnCounter; rowCounter++)
		{
			combinedRaw((columnCounter+nMarkers1)*(columnCounter+nMarkers1+1ULL)/2ULL + (rowCounter + nMarkers1)) = raw2(columnCounter*(columnCounter+1ULL)/2ULL + rowCounter);
		}
	}
	Rcpp::Function c("c");

	Rcpp::S4 retVal("rf");
	Rcpp::S4 theta("rawSymmetricMatrix");
	theta.slot("data") = combinedRaw;
	theta.slot("markers") = c(markers1, markers2);
	theta.slot("levels") = levels1;

	retVal.slot("theta") = theta;
	if(!Rcpp::as<Rcpp::RObject>(rf1.slot("lod")).isNULL() ||!Rcpp::as<Rcpp::RObject>(rf2.slot("lod")).isNULL())
	{
		throw std::runtime_error("Code for slot lod not implemented yet");
	}
	if(!Rcpp::as<Rcpp::RObject>(rf1.slot("lkhd")).isNULL() ||!Rcpp::as<Rcpp::RObject>(rf2.slot("lkhd")).isNULL())
	{
		throw std::runtime_error("Code for slot lkhd not implemented yet");
	}
	return retVal;
END_RCPP
}
