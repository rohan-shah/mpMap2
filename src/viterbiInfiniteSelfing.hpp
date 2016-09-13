#ifndef VITERBI_INFINITE_SELFING_HEADER_GUARD
#define VITERBI_INFINITE_SELFING_HEADER_GUARD
#include "intercrossingAndSelfingGenerations.h"
#include "recodeFoundersFinalsHets.h"
#include "matrices.hpp"
#include "probabilities.hpp"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#include "funnelsToUniqueValues.h"
#include "estimateRFCheckFunnels.h"
#include "markerPatternsToUniqueValues.h"
#include "intercrossingHaplotypeToMarker.hpp"
#include "funnelHaplotypeToMarker.hpp"
#include <limits>
template<int nFounders> struct viterbiAlgorithm<nFounders, true>
{
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	rowMajorMatrix<int> intermediate1, intermediate2;
	Rcpp::IntegerMatrix results;
	std::vector<double> pathLengths1, pathLengths2;
	std::vector<double> working;
	xMajorMatrix<expandedProbabilitiesType>* logIntercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>* logFunnelHaplotypeProbabilities;
	markerPatternsToUniqueValuesArgs& markerData;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minSelfingGenerations;
	int maxSelfingGenerations;
	int minAIGenerations, maxAIGenerations;
	double heterozygoteMissingProb, homozygoteMissingProb, errorProb;
	Rcpp::IntegerMatrix key;
	std::vector<array2<nFounders> >* logIntercrossingSingleLociHaplotypeProbabilities, *intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* logFunnelSingleLociHaplotypeProbabilities, *funnelSingleLociHaplotypeProbabilities;
	viterbiAlgorithm(markerPatternsToUniqueValuesArgs& markerData, int maxChromosomeSize)
		: intermediate1(nFounders, maxChromosomeSize), intermediate2(nFounders, maxChromosomeSize), pathLengths1(nFounders), pathLengths2(nFounders), working(nFounders), logIntercrossingHaplotypeProbabilities(NULL), logFunnelHaplotypeProbabilities(NULL), markerData(markerData), lineFunnelIDs(NULL), lineFunnelEncodings(NULL), intercrossingGenerations(NULL), selfingGenerations(NULL), minSelfingGenerations(-1), maxSelfingGenerations(-1), minAIGenerations(-1), maxAIGenerations(-1), heterozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), homozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), errorProb(std::numeric_limits<double>::quiet_NaN()), logIntercrossingSingleLociHaplotypeProbabilities(NULL), intercrossingSingleLociHaplotypeProbabilities(NULL), logFunnelSingleLociHaplotypeProbabilities(NULL), funnelSingleLociHaplotypeProbabilities(NULL)
	{}
	void apply(int start, int end)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || lineFunnelIDs == NULL || lineFunnelEncodings == NULL || intercrossingGenerations == NULL || selfingGenerations == NULL)
		{
			throw std::runtime_error("Internal error");
		}
		minAIGenerations = *std::min_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		maxAIGenerations = *std::max_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		minAIGenerations = std::max(minAIGenerations, 1);
		int nFinals = recodedFinals.nrow();
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if((*intercrossingGenerations)[finalCounter] == 0)
			{
				applyFunnel(start, end, finalCounter, (*lineFunnelIDs)[finalCounter]);
			}
			else
			{
				applyIntercrossing(start, end, finalCounter, (*intercrossingGenerations)[finalCounter]);
			}
			std::vector<double>::iterator longestPath = std::max_element(pathLengths1.begin(), pathLengths1.end());
			int longestIndex = (int)std::distance(pathLengths1.begin(), longestPath);
			for(int i = 0; i < end - start; i++)
			{
				results(finalCounter, i+start) = intermediate1(longestIndex, i);
			}
		}
	}
	void applyFunnel(int start, int end, int finalCounter, int funnelID)
	{
		if(errorProb == 0)
		{
			applyFunnelNoError(start, end, finalCounter, funnelID);
		}
		else
		{
			applyFunnelWithError(start, end, finalCounter, funnelID);
		}
	}
	void applyFunnelNoError(int start, int end, int finalCounter, int funnelID)
	{
		if(logFunnelHaplotypeProbabilities == NULL || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		int markerValue = recodedFinals(finalCounter, start);
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			intermediate1(founderCounter, 0) = intermediate2(founderCounter, 0) = founderCounter+1;
			pathLengths2[founderCounter] = pathLengths1[founderCounter] = 0;
			if(recodedFounders(founderCounter, start) != markerValue && markerValue != NA_INTEGER)
			{
				pathLengths2[founderCounter] = pathLengths1[founderCounter] = -std::numeric_limits<double>::infinity();
			}
		}
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(funnel[founderCounter], markerCounter+1) == markerValue || markerValue == NA_INTEGER)
				{
					//Founder at the previous marker. 
					std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						if(recodedFounders(funnel[founderCounter2], markerCounter) == previousMarkerValue || previousMarkerValue == NA_INTEGER)
						{
							working[funnel[founderCounter2]] = pathLengths1[funnel[founderCounter2]] + (*logFunnelHaplotypeProbabilities)(markerCounter-start, 0).values[founderCounter2][founderCounter];
						}
					}
					//Get the shortest one, and check that it's not negative infinity.
					std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
					int bestPrevious = (int)std::distance(working.begin(), longest);
					
					memcpy(&(intermediate2(funnel[founderCounter], identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(markerCounter - start + 1 - identicalIndex));
					intermediate2(funnel[founderCounter], markerCounter-start+1) = funnel[founderCounter]+1;
					pathLengths2[funnel[founderCounter]] = *longest;
				}
				else
				{
					pathLengths2[funnel[founderCounter]] = -std::numeric_limits<double>::infinity();
				}
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(0, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					if(value != intermediate1(founderCounter, identicalIndex)) goto stopIdenticalSearch;
				}
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					intermediate2(founderCounter, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyFunnelWithError(int start, int end, int finalCounter, int funnelID)
	{
		if(logFunnelHaplotypeProbabilities == NULL || errorProb == 0 || errorProb == 1)
		{
			throw std::runtime_error("Internal error");
		}
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		int markerValue = recodedFinals(finalCounter, start);
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
		//Table of the number of founders that map to that marker allele
		int table[nFounders];
		memset(table, 0, sizeof(int)*nFounders);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			table[startMarkerData.hetData(founderCounter, founderCounter)]++;
		}
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			intermediate1(founderCounter, 0) = intermediate2(founderCounter, 0) = founderCounter+1;
			pathLengths2[founderCounter] = pathLengths1[founderCounter] = 0;
			if(recodedFounders(founderCounter, start) == markerValue)
			{
				pathLengths2[founderCounter] = pathLengths1[founderCounter] = log((1.0 / (double)nFounders) * ((1 - errorProb) + errorProb * table[markerValue] / (double)nFounders));
			}
			else if(markerValue == NA_INTEGER)
			{
				pathLengths2[founderCounter] = pathLengths1[founderCounter] = log(1.0 / (double)nFounders);
			}
			else pathLengths2[founderCounter] = pathLengths1[founderCounter] = log((1.0 / (double)nFounders) * errorProb * table[markerValue] / (double)nFounders);
		}
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter]];
			memset(table, 0, sizeof(int)*nFounders);
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				table[currentMarkerData.hetData(founderCounter, founderCounter)]++;
			}
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				double increment;
				if(recodedFounders(funnel[founderCounter], markerCounter+1) == markerValue)
				{
					increment = log((1 - errorProb) + errorProb * table[markerValue] / (double)nFounders);
				}
				else if(markerValue == NA_INTEGER)
				{
					increment = 0;
				}
				else increment = log(errorProb * table[markerValue] / (double)nFounders);
				//Founder at the previous marker. 
				std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
				for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
				{
					working[funnel[founderCounter2]] = pathLengths1[funnel[founderCounter2]] + (*logFunnelHaplotypeProbabilities)(markerCounter-start, 0).values[founderCounter2][founderCounter] + increment;
				}
				//Get the shortest one, and check that it's not negative infinity.
				std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
				int bestPrevious = (int)std::distance(working.begin(), longest);
				
				memcpy(&(intermediate2(funnel[founderCounter], identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(markerCounter - start + 1 - identicalIndex));
				intermediate2(funnel[founderCounter], markerCounter-start+1) = funnel[founderCounter]+1;
				pathLengths2[funnel[founderCounter]] = *longest;
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(0, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					if(value != intermediate1(founderCounter, identicalIndex)) goto stopIdenticalSearch;
				}
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					intermediate2(founderCounter, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGeneration)
	{
		if(errorProb == 0)
		{
			applyIntercrossingNoError(start, end, finalCounter, intercrossingGeneration);
		}
		else
		{
			applyIntercrossingWithError(start, end, finalCounter, intercrossingGeneration);
		}
	}
	void applyIntercrossingNoError(int start, int end, int finalCounter, int intercrossingGeneration)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		int markerValue = recodedFinals(finalCounter, start);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			intermediate1(founderCounter, 0) = founderCounter+1;
			pathLengths1[founderCounter] = 0;
			if(recodedFounders(founderCounter, start) != markerValue && markerValue != NA_INTEGER)
			{
				pathLengths1[founderCounter] = -std::numeric_limits<double>::infinity();
			}
		}
		//The index, before which all the paths are identical
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				//NA corresponds to no restriction from the marker value
				if(recodedFounders(founderCounter, markerCounter+1) == markerValue || markerValue == NA_INTEGER)
				{
					//Founder at the previous marker. 
					std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						//NA corresponds to no restriction
						if(recodedFounders(founderCounter2, markerCounter) == previousMarkerValue || previousMarkerValue == NA_INTEGER)
						{
							working[founderCounter2] = pathLengths1[founderCounter2] + (*logIntercrossingHaplotypeProbabilities)(markerCounter-start, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
						}
					}
					//Get the longest one, and check that it's not negative infinity.
					std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
					int bestPrevious = (int)std::distance(working.begin(), longest);
					
					memcpy(&(intermediate2(founderCounter, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(markerCounter - start + 1 - identicalIndex));
					intermediate2(founderCounter, markerCounter-start+1) = founderCounter+1;
					pathLengths2[founderCounter] = *longest;
				}
				else
				{
					pathLengths2[founderCounter] = -std::numeric_limits<double>::infinity();
				}
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(0, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					if(value != intermediate1(founderCounter, identicalIndex)) goto stopIdenticalSearch;
				}
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					intermediate2(founderCounter, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossingWithError(int start, int end, int finalCounter, int intercrossingGeneration)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || errorProb == 0 || errorProb == 1)
		{
			throw std::runtime_error("Internal error");
		}
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
		//Table of the number of founders that map to that marker allele
		int table[nFounders];
		memset(table, 0, sizeof(int)*nFounders);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			table[startMarkerData.hetData(founderCounter, founderCounter)]++;
		}
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		int markerValue = recodedFinals(finalCounter, start);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			intermediate1(founderCounter, 0) = intermediate2(founderCounter, 0) = founderCounter+1;
			pathLengths2[founderCounter] = pathLengths1[founderCounter] = 0;
			if(recodedFounders(founderCounter, start) == markerValue)
			{
				pathLengths2[founderCounter] = pathLengths1[founderCounter] = log((1.0 / (double)nFounders) * ((1 - errorProb) + errorProb * table[markerValue] / (double)nFounders));
			}
			else if(markerValue == NA_INTEGER)
			{
				pathLengths2[founderCounter] = pathLengths1[founderCounter] = log(1.0 / (double)nFounders);
			}
			else pathLengths2[founderCounter] = pathLengths1[founderCounter] = log((1.0 / (double)nFounders) * errorProb * table[markerValue] / (double)nFounders);
		}
		//The index, before which all the paths are identical
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter]];
			memset(table, 0, sizeof(int)*nFounders);
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				table[currentMarkerData.hetData(founderCounter, founderCounter)]++;
			}
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				//NA corresponds to no restriction from the marker value
				double increment;
				if(recodedFounders(founderCounter, markerCounter+1) == markerValue)
				{
					increment = log((1 - errorProb) + errorProb * table[markerValue] / (double)nFounders);;
				}
				else if(markerValue == NA_INTEGER)
				{
					increment = 0;
				}
				else increment = log(errorProb * table[markerValue] / (double)nFounders);
				//Founder at the previous marker. 
				std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
				for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
				{
					working[founderCounter2] = pathLengths1[founderCounter2] + (*logIntercrossingHaplotypeProbabilities)(markerCounter-start, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter] + increment;
				}
				//Get the longest one, and check that it's not negative infinity.
				std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
				int bestPrevious = (int)std::distance(working.begin(), longest);
				
				memcpy(&(intermediate2(founderCounter, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(markerCounter - start + 1 - identicalIndex));
				intermediate2(founderCounter, markerCounter-start+1) = founderCounter+1;
				pathLengths2[founderCounter] = *longest;
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(0, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					if(value != intermediate1(founderCounter, identicalIndex)) goto stopIdenticalSearch;
				}
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					intermediate2(founderCounter, identicalIndex) = value;
				}
				identicalIndex++;
			}
	stopIdenticalSearch:
			;
		}
	}
};
#endif
