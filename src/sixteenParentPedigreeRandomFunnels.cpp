#include "sixteenParentPedigreeRandomFunnels.h"
#include "Rcpp.h"
//pairs <- combn(1:16, 2)
//pairs <- apply(pairs, 2, function(x) c(min(x), max(x)))
//apply(pairs, 2, function(x) cat("{", x[1], ", ", x[2], "}, ", sep=""))
//120 rows, 2 column
const int sixteenParentRandomFunnelPairs[][2] = {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}, {1, 11}, {1, 12}, {1, 13}, {1, 14}, {1, 15}, {1, 16}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {2, 7}, {2, 8}, {2, 9}, {2, 10}, {2, 11}, {2, 12}, {2, 13}, {2, 14}, {2, 15}, {2, 16}, {3, 4}, {3, 5}, {3, 6}, {3, 7}, {3, 8}, {3, 9}, {3, 10}, {3, 11}, {3, 12}, {3, 13}, {3, 14}, {3, 15}, {3, 16}, {4, 5}, {4, 6}, {4, 7}, {4, 8}, {4, 9}, {4, 10}, {4, 11}, {4, 12}, {4, 13}, {4, 14}, {4, 15}, {4, 16}, {5, 6}, {5, 7}, {5, 8}, {5, 9}, {5, 10}, {5, 11}, {5, 12}, {5, 13}, {5, 14}, {5, 15}, {5, 16}, {6, 7}, {6, 8}, {6, 9}, {6, 10}, {6, 11}, {6, 12}, {6, 13}, {6, 14}, {6, 15}, {6, 16}, {7, 8}, {7, 9}, {7, 10}, {7, 11}, {7, 12}, {7, 13}, {7, 14}, {7, 15}, {7, 16}, {8, 9}, {8, 10}, {8, 11}, {8, 12}, {8, 13}, {8, 14}, {8, 15}, {8, 16}, {9, 10}, {9, 11}, {9, 12}, {9, 13}, {9, 14}, {9, 15}, {9, 16}, {10, 11}, {10, 12}, {10, 13}, {10, 14}, {10, 15}, {10, 16}, {11, 12}, {11, 13}, {11, 14}, {11, 15}, {11, 16}, {12, 13}, {12, 14}, {12, 15}, {12, 16}, {13, 14}, {13, 15}, {13, 16}, {14, 15}, {14, 16}, {15, 16}};
//pairIndices <- matrix(NA, 16, 16)
//for(i in 1:16) for(j in setdiff(1:16, i)) pairIndices[i, j] <- which(pairs[1,] == min(i, j) & pairs[2,] == max(i, j))
//pairIndices[is.na(pairIndices)] <- -1000
//apply(pairIndices, 1, function(x) cat("{", do.call(paste, c(as.list(x), sep=", ")), "}, ", sep=""))
const int sixteenParentRandomFunnelPairIndices[][16] = {{-1000, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {1, -1000, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}, {2, 16, -1000, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42}, {3, 17, 30, -1000, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54}, {4, 18, 31, 43, -1000, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65}, {5, 19, 32, 44, 55, -1000, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75}, {6, 20, 33, 45, 56, 66, -1000, 76, 77, 78, 79, 80, 81, 82, 83, 84}, {7, 21, 34, 46, 57, 67, 76, -1000, 85, 86, 87, 88, 89, 90, 91, 92}, {8, 22, 35, 47, 58, 68, 77, 85, -1000, 93, 94, 95, 96, 97, 98, 99}, {9, 23, 36, 48, 59, 69, 78, 86, 93, -1000, 100, 101, 102, 103, 104, 105}, {10, 24, 37, 49, 60, 70, 79, 87, 94, 100, -1000, 106, 107, 108, 109, 110}, {11, 25, 38, 50, 61, 71, 80, 88, 95, 101, 106, -1000, 111, 112, 113, 114}, {12, 26, 39, 51, 62, 72, 81, 89, 96, 102, 107, 111, -1000, 115, 116, 117}, {13, 27, 40, 52, 63, 73, 82, 90, 97, 103, 108, 112, 115, -1000, 118, 119}, {14, 28, 41, 53, 64, 74, 83, 91, 98, 104, 109, 113, 116, 118, -1000, 120}, {15, 29, 42, 54, 65, 75, 84, 92, 99, 105, 110, 114, 117, 119, 120, -1000}};

SEXP sixteenParentPedigreeRandomFunnels(SEXP initialPopulationSize_sexp, SEXP selfingGenerations_sexp, SEXP nSeeds_sexp, SEXP intercrossingGenerations_sexp)
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
	Rcpp::IntegerMatrix funnels(initialPopulationSize, 16);
	Rcpp::Function sample("sample"), setdiff("setdiff");
	for(int i = 0; i < initialPopulationSize; i++)
	{
		Rcpp::IntegerVector sampled = sample(Rcpp::Range(1, 16), 16, Rcpp::Named("replace") = false);
		for(int j = 0; j < 16; j++) funnels(i, j) = sampled(j);
	}
	int entries = 16 + 120 + 7 * initialPopulationSize + intercrossingGenerations * initialPopulationSize + nSeeds * selfingGenerations * initialPopulationSize;
	Rcpp::IntegerVector mother(entries), father(entries);
	Rcpp::LogicalVector observed(entries, false);

	for(int i = 0; i < 16; i++)
	{
		mother(i) = father(i) = 0;
	}

	//Put in all pairs
	for(int i = 0; i < 120; i++)
	{
		mother[16 + i] = sixteenParentRandomFunnelPairs[i][0];
		father[16 + i] = sixteenParentRandomFunnelPairs[i][1];
		if(sixteenParentRandomFunnelPairs[i][0] < 0 || sixteenParentRandomFunnelPairs[i][1] < 0)
		{
			throw std::runtime_error("Internal error");
		}
	}

	//We have to generate new funnels for every line
	for(int i = 0; i < initialPopulationSize; i++)
	{
		for(int pairIndex = 0; pairIndex < 4; pairIndex++)
		{
			mother[16 + 120 + i*4 + pairIndex] = sixteenParentRandomFunnelPairIndices[funnels(i, 4*pairIndex)-1][funnels(i, 4*pairIndex+1)-1]+16;
			father[16 + 120 + i*4 + pairIndex] = sixteenParentRandomFunnelPairIndices[funnels(i, 4*pairIndex+2)-1][funnels(i, 4*pairIndex+3)-1]+16;
		}
	}
	//Create two eight-way crosses for each funnel
	for(int i = 0; i < 2*initialPopulationSize; i++)
	{
		mother[16 + 120 + 4*initialPopulationSize + i] = 15 + 120 + 2 * (i+1);
		father[16 + 120 + 4*initialPopulationSize + i] = 16 + 120 + 2 * (i+1);
	}
	
	//Create a line which mixes all sixteen founders, for each funnel
	for(int i = 0; i < initialPopulationSize; i++)
	{
		mother[16 + 120 + 6*initialPopulationSize + i] = 15 + 120 + 4*initialPopulationSize + 2 * (i+1);
		father[16 + 120 + 6*initialPopulationSize + i] = 16 + 120 + 4*initialPopulationSize + 2 * (i+1);
	}

	int currentIndex = 16 + 120 + initialPopulationSize*6;
	int lastGenerationStart = currentIndex;
	int lastGenerationEnd = currentIndex + initialPopulationSize;
	if(intercrossingGenerations > 0)
	{
		for(int intecrossingGenerationCounter = 0; intecrossingGenerationCounter < intercrossingGenerations; intecrossingGenerationCounter++)
		{
			for(int lineCounter = lastGenerationStart; lineCounter < lastGenerationEnd; lineCounter++)
			{
				mother(lineCounter + initialPopulationSize) = lineCounter+1;
				Rcpp::IntegerVector possibilities = setdiff(Rcpp::Range(lastGenerationStart, lastGenerationEnd-1), lineCounter);
				father(lineCounter + initialPopulationSize) = Rcpp::as<int>(sample(possibilities, 1))+1;
			}
			lastGenerationStart += initialPopulationSize;
			lastGenerationEnd += initialPopulationSize;
		}
		currentIndex = lastGenerationStart;
	}
	int nextFree = currentIndex + initialPopulationSize;
	if(selfingGenerations == 1)
	{
		for(int lineCounter = currentIndex; lineCounter < currentIndex + initialPopulationSize; lineCounter++)
		{
			for(int seedCounter = nextFree; seedCounter < nextFree + nSeeds; seedCounter++)
			{
				mother(seedCounter) = father(seedCounter) = lineCounter+1;
			}
			observed(nextFree+nSeeds-1) = true;
			nextFree = nextFree + nSeeds;
		}
	}
	else if(selfingGenerations > 1)
	{
		for(int lineCounter = currentIndex; lineCounter < currentIndex + initialPopulationSize; lineCounter++)
		{
			for(int seedCounter = 0; seedCounter < nSeeds; seedCounter++)
			{
				father(nextFree) = mother(nextFree) = lineCounter+1;
				for(int j = nextFree+1; j < nextFree + selfingGenerations; j++)
				{
					father(j) = mother(j) = j;
				}
				observed(nextFree+selfingGenerations-1) = true;
				nextFree += selfingGenerations;
			}
		}
	}
	else
	{
		if(intercrossingGenerations == 0)
		{
			for(int i = nextFree-initialPopulationSize+1; i < entries; i++)
			{
				observed[i] = true;
			}
		}
		else
		{
			for(int i = lastGenerationStart; i < lastGenerationEnd; i++)
			{
				observed[i] = true;
			}
		}
	}
	Rcpp::Function newCall("new");
	Rcpp::Function paste0("paste0");
	Rcpp::CharacterVector lineNames = paste0("L", Rcpp::Range(1, entries));

	Rcpp::S4 result = newCall("detailedPedigree", Rcpp::Named("lineNames") = lineNames, Rcpp::Named("mother") = mother, Rcpp::Named("father") = father, Rcpp::Named("initial") = Rcpp::Range(1, 16), Rcpp::Named("observed") = observed, Rcpp::Named("selfing") = "infinite", Rcpp::Named("warnImproperFunnels") = true);
	return result;
END_RCPP
}
