#include "fourParentPedigreeSingleFunnel.h"
SEXP fourParentPedigreeSingleFunnel(SEXP initialPopulationSize_sexp, SEXP selfingGenerations_sexp, SEXP nSeeds_sexp, SEXP intercrossingGenerations_sexp)
{
BEGIN_RCPP
	int initialPopulationSize, selfingGenerations, nSeeds, intercrossingGenerations;
	try
	{
		initialPopulationSize = Rcpp::as<int>(initialPopulationSize_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Argument initialPopulationSize must be an integer");
	}
	try
	{
		selfingGenerations = Rcpp::as<int>(selfingGenerations_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Argument selfingGenerations must be an integer");
	}
	try
	{
		nSeeds = Rcpp::as<int>(nSeeds_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Argument nSeeds must be an integer");
	}
	try
	{
		intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Argument intercrossingGenerations must be an integer");
	}
	int nEntries = 4 + 2 + initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize;

	//R functions that we're going to call
	Rcpp::Function paste0("paste0"), sample("sample"), setdiff("setdiff"), newCall("new");
	
	Rcpp::IntegerVector mother(nEntries, NA_INTEGER), father(nEntries, NA_INTEGER);
	Rcpp::LogicalVector observed(nEntries, false);
	Rcpp::CharacterVector lineNames = paste0("L", Rcpp::Range(1, nEntries));

	std::fill(mother.begin(), mother.begin()+4, 0);
	std::fill(father.begin(), father.begin()+4, 0);

	mother(4) = 1;
	mother(5) = 3;

	father(4) = 2;
	father(5) = 4;

	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother(i+6) = 5;
		father(i+6) = 6;
	}

	int currentIndex = 6;
	int lastGenerationStart = currentIndex;
	int lastGenerationEnd = currentIndex + initialPopulationSize;
	if(intercrossingGenerations > 0)
	{
		for(int i = 0; i < intercrossingGenerations; i++)
		{
			Rcpp::IntegerVector possibilities = Rcpp::Range(lastGenerationStart+1, lastGenerationEnd-1);
			for(int lineCounter = lastGenerationStart; lineCounter < lastGenerationEnd; lineCounter++)
			{
				mother(lineCounter + initialPopulationSize) = lineCounter+1;

				Rcpp::IntegerVector sampled = sample(possibilities, 1);
				father(lineCounter + initialPopulationSize) = Rcpp::as<int>(sampled) + 1;
				if(lineCounter != lastGenerationEnd - 1) possibilities(lineCounter - lastGenerationStart) = lineCounter;
			}
			lastGenerationStart += initialPopulationSize;
			lastGenerationEnd += initialPopulationSize;
		}
		currentIndex = lastGenerationStart;
	}
	int nextFree = currentIndex+initialPopulationSize;
	if(selfingGenerations == 1)
	{
		for(int lineCounter = currentIndex; lineCounter < currentIndex+initialPopulationSize; lineCounter++)
		{
			std::fill(mother.begin() + nextFree, mother.begin() + nextFree + nSeeds, lineCounter + 1);
			std::fill(father.begin() + nextFree, father.begin() + nextFree + nSeeds, lineCounter + 1);
			observed(nextFree+nSeeds-1) = true;
			nextFree += nSeeds;
		}
	}
	else if(selfingGenerations > 1)
	{
		for(int lineCounter = currentIndex; lineCounter < currentIndex+initialPopulationSize; lineCounter++)
		{
			for(int seedCounter = 0; seedCounter < nSeeds; seedCounter++)
			{
				father(nextFree) = mother(nextFree) = lineCounter + 1;
				for(int i = nextFree + 1; i < nextFree+selfingGenerations; i++)
				{
					father(i) = mother(i) = i;
				}
				observed(nextFree + selfingGenerations-1) = true;
				nextFree += selfingGenerations;
			}
		}
	}
	else
	{
		if(intercrossingGenerations == 0)
		{
			std::fill(observed.begin()+6, observed.end(), true);
		}
		else
		{
			std::fill(observed.begin()+lastGenerationStart, observed.begin() + lastGenerationEnd, true);
		}
	}
	Rcpp::S4 result = newCall("detailedPedigree", Rcpp::Named("lineNames") = lineNames, Rcpp::Named("mother") = mother, Rcpp::Named("father") = father, Rcpp::Named("initial") = Rcpp::Range(1, 4), Rcpp::Named("observed") = observed, Rcpp::Named("selfing") = "infinite", Rcpp::Named("warnImproperFunnels") = true);
	return result;
END_RCPP
}
