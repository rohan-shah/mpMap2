#include "generateGenotypes.h"
#include <Rcpp.h>
void createGamete(Rcpp::NumericVector& recombinationFractions, Rcpp::IntegerVector& geneticData, int* output)
{
	R_xlen_t nMarkers = recombinationFractions.length() + 1;
	//Set the first marker separately as it's set without reference to any other
	int runningHaplotype = 0;
	if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) < 0.5) runningHaplotype = 1;
	output[0] = geneticData[nMarkers * runningHaplotype];
	GetRNGstate();
	for(int i = 1; i < nMarkers; i++)
	{
		//recombination if we fall below the relevant recombination fraction
		if(::unif_rand() < recombinationFractions(i-1)) runningHaplotype = abs(runningHaplotype-1);
		output[i] = geneticData[nMarkers * runningHaplotype + i];
	}
	PutRNGstate();
}
SEXP generateGenotypes(SEXP RrecombinationFractions, SEXP RmarkerNames, SEXP Rpedigree)
{
	BEGIN_RCPP
		Rcpp::NumericVector recombinationFractions;
		try
		{
			recombinationFractions = Rcpp::NumericVector(RrecombinationFractions);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input recombinationFractions must be a numeric vector");
		}
		for(std::size_t i = 0; i < (std::size_t)recombinationFractions.size(); i++)
		{
			if(recombinationFractions[i] < 0 || recombinationFractions[i] > 0.5) throw std::runtime_error("All recombination fractions must be between 0 and 0.5, inclusive");
		}
		Rcpp::S4 pedigree;
		try
		{
			pedigree = Rcpp::S4(Rpedigree);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input pedigree must be an S4 object");
		}
		Rcpp::CharacterVector markerNames;
		try
		{
			markerNames = RmarkerNames;
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input markerNames must be a character vector");
		}
		if(Rcpp::as<std::string>(pedigree.attr("class")) != "detailedPedigree")
		{
			throw std::runtime_error("Input pedigree had class " + Rcpp::as<std::string>(pedigree.attr("class")) + " instead of detailedPedigree");
		}
		Rcpp::CharacterVector lineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));
		Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother"));
		Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));
		
		R_xlen_t nMarkers = recombinationFractions.length() + 1, nPedRows = lineNames.size();
		if(nMarkers != markerNames.size()) throw std::runtime_error("Input marker names vector had the wrong length");
		//Columns 0 and nMarkers correspond to the the pair of alleles at the first marker
		Rcpp::IntegerMatrix result((int)lineNames.size(), (int)(2*nMarkers));
		//once we encounter a line with parents which are not set to zero we set this to true. And subsequently don't allow any zero-values for mother or father
		bool finishedFounders = false;
		int founderCounter = 0;
		//count number of founders by looking at the number of rows with mother and founder both set to '0'. This DOESN'T have to be a power of 2. It just indicates the number of lines which are generated without reference to any parent lines.
		//So they're generated a little differently. 
		for(int lineCounter = 0; lineCounter < nPedRows; lineCounter++)
		{
			bool isFounder = (mother(lineCounter) == 0 || father(lineCounter) == 0);
			if(finishedFounders && (isFounder))
			{
				throw std::runtime_error("Founder lines must be placed at the top of the pedigree");
			}
			if(isFounder)
			{
				//If it's a founder, fill it with the index of the founder
				for(int counter = 0; counter < 2*nMarkers; counter++)
				{
					result(lineCounter, counter) = founderCounter+1;
				}
				founderCounter++;
			}
			else
			{
				//otherwise use the previously generated genetic data for the mother and father
				int motherPedigreeRowIndex = mother(lineCounter);
				int fatherPedigreeRowIndex = father(lineCounter);
				//copy out the genetic data for mother and father (includes BOTH alleles at every location)
				Rcpp::IntegerVector motherData = result.row(motherPedigreeRowIndex-1);
				Rcpp::IntegerVector fatherData = result.row(fatherPedigreeRowIndex-1);
				
				Rcpp::IntegerVector newGenotypes(2*nMarkers);
				createGamete(recombinationFractions, motherData, &(newGenotypes(0)));
				createGamete(recombinationFractions, fatherData, &(newGenotypes(0)) + nMarkers);
				result.row(lineCounter) = newGenotypes;

				finishedFounders = true;
			}
		}
		Rcpp::CharacterVector resultColNames(2*nMarkers);
		for(int i = 0; i < nMarkers; i++) resultColNames[i] = resultColNames[i+nMarkers] = markerNames[i];
		Rcpp::List resultDimNames = Rcpp::List::create(lineNames, resultColNames);

		result.attr("dimnames") = resultDimNames;
		return result;
	END_RCPP
}
