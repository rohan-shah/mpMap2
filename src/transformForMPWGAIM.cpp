#include "transformForMPWGAIM.h"
#include "recodeFoundersFinalsHets.h"
SEXP transformForMPWGAIM(SEXP geneticData_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData = Rcpp::as<Rcpp::S4>(geneticData_sexp);
	Rcpp::IntegerMatrix finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	Rcpp::IntegerMatrix founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	Rcpp::S4 hetData = Rcpp::as<Rcpp::S4>(geneticData.slot("hetData"));
	Rcpp::S4 probabilities = Rcpp::as<Rcpp::S4>(geneticData.slot("probabilities"));
	Rcpp::NumericMatrix probabilitiesData = Rcpp::as<Rcpp::NumericMatrix>(probabilities.slot("data"));

	int nFounders = (int)founders.nrow();
	int nMarkers = (int)founders.ncol();
	int nFinals = (int)finals.nrow();

	if(probabilitiesData.nrow() != nFinals * nFounders || probabilitiesData.ncol() != nMarkers)
	{
		throw std::runtime_error("Probabilities data had the wrong dimensions");
	}
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	Rcpp::List recodedHetData(nMarkers);
	Rcpp::CharacterVector markerNames = hetData.attr("names");
	Rcpp::CharacterVector founderNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(founders.attr("dimnames"))(0));
	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))(0));

	recodedHetData.attr("names") = markerNames;
	recodedFinals.attr("dimnames") = finals.attr("dimnames");
	recodedFounders.attr("dimnames") = founders.attr("dimnames");

	recodeDataStruct recoded;
	recoded.recodedFounders = recodedFounders;
	recoded.recodedFinals = recodedFinals;
	recoded.recodedHetData = recodedHetData;
	recoded.founders = founders;
	recoded.finals = finals;
	recoded.hetData = hetData;
	recodeFoundersFinalsHets(recoded);

	bool hasHets = false, hasHetEncodings = false;
	for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		int maxEncoding = recodedFounders(0, markerCounter);
		for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
		{
			maxEncoding = std::max(maxEncoding, recodedFounders(founderCounter, markerCounter));
		}
		Rcpp::IntegerMatrix currentMarkerHetData = Rcpp::as<Rcpp::IntegerMatrix>(recodedHetData(markerCounter));
		for(int i = 0; i < currentMarkerHetData.nrow(); i++)
		{
			if(currentMarkerHetData(i, 2) > maxEncoding) hasHetEncodings = true;
		}
		bool isBiallelic = maxEncoding == 1;
		//In the biallelic case code the markers as -1 / +1
		if(isBiallelic)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				recodedFounders(founderCounter, markerCounter) = 2*recodedFounders(founderCounter, markerCounter)-1;
			}
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				if(recodedFinals(finalCounter, markerCounter) == NA_INTEGER)
				{}
				else if(recodedFinals(finalCounter, markerCounter) > maxEncoding)
				{
					recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
					hasHets = true;
				}
				else
				{
					recodedFinals(finalCounter, markerCounter) = 2*recodedFinals(finalCounter, markerCounter)-1;
				}
			}
		}
		//Otherwise keep the previous encoding
		else
		{
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				if(recodedFinals(finalCounter, markerCounter) == NA_INTEGER)
				{}
				else if(recodedFinals(finalCounter, markerCounter) > maxEncoding)
				{
					recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
					hasHets = true;
				}
			}
		}
	}
	//Probability data also has to be transformed
	Rcpp::NumericMatrix transformedProbabilities(nFinals, nMarkers*nFounders);
	Rcpp::CharacterVector transformedProbabilitiesColNames(nMarkers*nFounders);
	int counter = 0;
	for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			transformedProbabilitiesColNames(counter) = Rcpp::as<std::string>(markerNames(markerCounter)) + " - " + Rcpp::as<std::string>(founderNames(founderCounter));
			counter++;
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				transformedProbabilities(finalCounter, markerCounter*nFounders + founderCounter) = probabilitiesData(finalCounter * nFounders + founderCounter, markerCounter);
			}
		}
	}
	transformedProbabilities.attr("dimnames") = Rcpp::List::create(finalNames, transformedProbabilitiesColNames);
	Rcpp::List retVal = Rcpp::List::create(Rcpp::Named("finals") = recodedFinals, Rcpp::Named("founders") = recodedFounders, Rcpp::Named("probabilities") = transformedProbabilities);
	return retVal;
END_RCPP
}
