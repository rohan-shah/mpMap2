#include "eightParentPedigreeRandomFunnels.h"
SEXP eightParentPedigreeRandomFunnels(SEXP initialPopulationSize_sexp, SEXP selfingGenerations_sexp, SEXP nSeeds_sexp, SEXP intercrossingGenerations_sexp)
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
	int entries = 8 + 7 * initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize;

	//The funnels matrix could be hard-coded, but I'll stick to doing it this way for now.
	Rcpp::IntegerMatrix pairs(2, 28);
	int counter = 0;
	for(int i = 1; i < 9; i++)
	{
		for(int j = i+1; j < 9; j++)
		{
			pairs(0, counter) = i;
			pairs(1, counter) = j;
			counter++;
		}
	}

	Rcpp::IntegerMatrix quads(4, 210);
	counter = 0;
	std::vector<int> possibleQuad(4), copiedPossibleQuad(4);
	for(int pairCounter1 = 0; pairCounter1 < 28; pairCounter1++)
	{
		possibleQuad[0] = pairs(0, pairCounter1);
		possibleQuad[1] = pairs(1, pairCounter1);
		for(int pairCounter2 = pairCounter1 + 1; pairCounter2 < 28; pairCounter2++)
		{
			possibleQuad[2] = pairs(0, pairCounter2);
			possibleQuad[3] = pairs(1, pairCounter2);
			copiedPossibleQuad = possibleQuad;
			std::sort(copiedPossibleQuad.begin(), copiedPossibleQuad.end());
			if(std::unique(copiedPossibleQuad.begin(), copiedPossibleQuad.end()) == copiedPossibleQuad.end())
			{
				quads(0, counter) = possibleQuad[0];
				quads(1, counter) = possibleQuad[1];
				quads(2, counter) = possibleQuad[2];
				quads(3, counter) = possibleQuad[3];
				counter++;
			}
		}
	}

	Rcpp::IntegerMatrix funnels(8, 315);
	counter = 0;
	std::vector<int> possibleEight(8), copiedPossibleEight(8);
	for(int quadCounter1 = 0; quadCounter1 < 210; quadCounter1++)
	{
		possibleEight[0] = quads(0, quadCounter1);
		possibleEight[1] = quads(1, quadCounter1);
		possibleEight[2] = quads(2, quadCounter1);
		possibleEight[3] = quads(3, quadCounter1);
		for(int quadCounter2 = quadCounter1 + 1; quadCounter2 < 210; quadCounter2++)
		{
			possibleEight[4] = quads(0, quadCounter2);
			possibleEight[5] = quads(1, quadCounter2);
			possibleEight[6] = quads(2, quadCounter2);
			possibleEight[7] = quads(3, quadCounter2);
			copiedPossibleEight = possibleEight;
			std::sort(copiedPossibleEight.begin(), copiedPossibleEight.end());
			if(std::unique(copiedPossibleEight.begin(), copiedPossibleEight.end()) == copiedPossibleEight.end())
			{
				funnels(0, counter) = possibleEight[0];
				funnels(1, counter) = possibleEight[1];
				funnels(2, counter) = possibleEight[2];
				funnels(3, counter) = possibleEight[3];
				funnels(4, counter) = possibleEight[4];
				funnels(5, counter) = possibleEight[5];
				funnels(6, counter) = possibleEight[6];
				funnels(7, counter) = possibleEight[7];
				counter++;
			}
		}
	}
	Rcpp::IntegerVector mother(entries), father(entries);
	for(int i = 0; i < 8; i++) mother[i] = father[i] = 0;

	Rcpp::LogicalVector observed(entries, false);

	Rcpp::Function paste0("paste0");
	Rcpp::CharacterVector lineNames = paste0("L", Rcpp::Range(1, entries));

	Rcpp::Function sample("sample");
	Rcpp::IntegerVector funnelNumbers = sample(Rcpp::Range(0, 314), initialPopulationSize, Rcpp::Named("replace") = true);

	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother[8 + i*4 + 0] = funnels(0, funnelNumbers(i));
		mother[8 + i*4 + 1] = funnels(2, funnelNumbers(i));
		mother[8 + i*4 + 2] = funnels(4, funnelNumbers(i));
		mother[8 + i*4 + 3] = funnels(6, funnelNumbers(i));
		father[8 + i*4 + 0] = funnels(1, funnelNumbers(i));
		father[8 + i*4 + 1] = funnels(3, funnelNumbers(i));
		father[8 + i*4 + 2] = funnels(5, funnelNumbers(i));
		father[8 + i*4 + 3] = funnels(7, funnelNumbers(i));
	}

	for(int i = 0; i < 2*initialPopulationSize; i++)
	{
		mother(8 + initialPopulationSize*4 + i) = 7L + 2L*(i+1);
		father(8 + initialPopulationSize*4 + i) = 8L + 2L*(i+1);
	}
	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother(8 + initialPopulationSize*6 + i) = 7L + initialPopulationSize*4 + 2L*(i+1);
		father(8 + initialPopulationSize*6 + i) = 8L + initialPopulationSize*4 + 2L*(i+1);
	}

	int currentIndex = 8 + initialPopulationSize*5 + initialPopulationSize;
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
			std::fill(observed.begin()+nextFree-initialPopulationSize, observed.end(), true);
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
