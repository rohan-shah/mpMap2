#include "eightParentPedigreeSingleFunnel.h"
SEXP eightParentPedigreeSingleFunnel(SEXP initialPopulationSize_sexp, SEXP selfingGenerations_sexp, SEXP nSeeds_sexp, SEXP intercrossingGenerations_sexp)
{
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
	int entries = 12 + 3 * initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize;

	Rcpp::IntegerVector mother(entries), father(entries);
	for(int i = 0; i < 8; i++) mother[i] = father[i] = 0;

	mother(8) = 1;
	mother(9) = 3;
	mother(10) = 5;
	mother(11) = 7;

	father(8) = 2;
	father(9) = 4;
	father(10) = 6;
	father(11) = 8;

	Rcpp::LogicalVector observed(entries, false);

	Rcpp::Function paste0("paste0");
	Rcpp::CharacterVector lineNames = paste0("L", Rcpp::Range(1, entries));

	Rcpp::Function sample("sample");

	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother(12 + 2*i) = 9;
		mother(12 + 2*i + 1) = 11;
		father(12 + 2*i) = 10;
		father(12 + 2*i + 1) = 12;
	}
	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother(12 + initialPopulationSize*2 + i) = 13 + 2*i;
		father(12 + initialPopulationSize*2 + i) = 14 + 2*i;
	}

	int currentIndex = 12 + initialPopulationSize*2;
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
			std::fill(observed.begin()+nextFree-initialPopulationSize+1, observed.end(), true);
		}
		else
		{
			std::fill(observed.begin()+lastGenerationStart, observed.begin() + lastGenerationEnd, true);
		}
	}
	Rcpp::Function newCall("new");
	Rcpp::S4 result = newCall("detailedPedigree", Rcpp::Named("lineNames") = lineNames, Rcpp::Named("mother") = mother, Rcpp::Named("father") = father, Rcpp::Named("initial") = Rcpp::Range(1, 8), Rcpp::Named("observed") = observed, Rcpp::Named("selfing") = "infinite", Rcpp::Named("warnImproperFunnels") = true);
	return result;
}
