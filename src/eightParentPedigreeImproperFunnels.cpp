#include "eightParentPedigreeRandomFunnels.h"
SEXP eightParentPedigreeImproperFunnels(SEXP initialPopulationSize_sexp, SEXP selfingGenerations_sexp, SEXP nSeeds_sexp)
{
	int initialPopulationSize, selfingGenerations, nSeeds;
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
	int entries = 8 + 7 * initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize;

	Rcpp::Function sample("sample");
	Rcpp::IntegerVector funnels = sample(Rcpp::Range(1, 8), initialPopulationSize*8, Rcpp::Named("replace") = true);

	Rcpp::IntegerVector mother(entries), father(entries);
	for(int i = 0; i < 8; i++) mother[i] = father[i] = 0;

	Rcpp::LogicalVector observed(entries, false);

	Rcpp::Function paste0("paste0");
	Rcpp::CharacterVector lineNames = paste0("L", Rcpp::Range(1, entries));

	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother[8 + i*4 + 0] = funnels(i + initialPopulationSize*0);
		mother[8 + i*4 + 1] = funnels(i + initialPopulationSize*2);
		mother[8 + i*4 + 2] = funnels(i + initialPopulationSize*4);
		mother[8 + i*4 + 3] = funnels(i + initialPopulationSize*6);
		father[8 + i*4 + 0] = funnels(i + initialPopulationSize*1);
		father[8 + i*4 + 1] = funnels(i + initialPopulationSize*3);
		father[8 + i*4 + 2] = funnels(i + initialPopulationSize*5);
		father[8 + i*4 + 3] = funnels(i + initialPopulationSize*7);
	}

	for(int i = 0; i < 2*initialPopulationSize; i++)
	{
		mother(8 + initialPopulationSize*4 + i) = 7 + 2*(i+1);
		father(8 + initialPopulationSize*4 + i) = 8 + 2*(i+1);
	}
	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother(8 + initialPopulationSize*6 + i) = 7 + initialPopulationSize*4 + 2*(i+1);
		father(8 + initialPopulationSize*6 + i) = 8 + initialPopulationSize*4 + 2*(i+1);
	}

	int currentIndex = 8 + initialPopulationSize*5 + initialPopulationSize;
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
		std::fill(observed.begin()+nextFree-initialPopulationSize+1, observed.end(), true);
	}
	Rcpp::Function newCall("new");
	Rcpp::S4 result = newCall("detailedPedigree", Rcpp::Named("lineNames") = lineNames, Rcpp::Named("mother") = mother, Rcpp::Named("father") = father, Rcpp::Named("initial") = Rcpp::Range(1, 8), Rcpp::Named("observed") = observed, Rcpp::Named("selfing") = "infinite", Rcpp::Named("warnImproperFunnels") = false);
	return result;
}
