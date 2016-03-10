#ifndef VITERBI_FINITE_SELFING_HEADER_GUARD
#define VITERBI_FINITE_SELFING_HEADER_GUARD
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
template<int nFounders> struct viterbiAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	rowMajorMatrix<int> intermediate1, intermediate2;
	Rcpp::IntegerMatrix results;
	std::vector<double> pathLengths1, pathLengths2;
	std::vector<double> working;
	xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities;
	markerPatternsToUniqueValuesArgs& markerData;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minSelfingGenerations;
	int maxSelfingGenerations;
	int minAIGenerations;
	int maxAIGenerations;
	Rcpp::IntegerMatrix key;
	double hetrozygoteMissingProb, homozygoteMissingProb;
	std::vector<array2<nFounders> >* intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* funnelSingleLociHaplotypeProbabilities;
	viterbiAlgorithm(markerPatternsToUniqueValuesArgs& markerData, xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities, rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities, int maxChromosomeSize)
		: markerData(markerData), intermediate1(nFounders*nFounders, maxChromosomeSize), intermediate2(nFounders*nFounders, maxChromosomeSize), pathLengths1(nFounders*nFounders), pathLengths2(nFounders*nFounders), working(nFounders*nFounders), intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities)
	{}
	void apply(int start, int end)
	{
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
		minAIGenerations = *std::min_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		maxAIGenerations = *std::max_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		minAIGenerations = std::max(minAIGenerations, 1);
		int nFinals = recodedFinals.nrow(), nMarkers = recodedFinals.ncol();

		//If there's not meant to be any missing values, check that first
		if(homozygoteMissingProb == 0 && hetrozygoteMissingProb == 0)
		{
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				for(int markerCounter = start; markerCounter < end; markerCounter++)
				{
					if(recodedFinals(finalCounter, markerCounter) == NA_INTEGER)
					{
						throw std::runtime_error("Inputs hetrozygoteMissingProb and homozygoteMissingProb imply that missing values are not allowed");
					}
				}
			}
		}
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if((*intercrossingGenerations)[finalCounter] == 0)
			{
				applyFunnel(start, end, finalCounter, (*lineFunnelIDs)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			else
			{
				applyIntercrossing(start, end, finalCounter, (*intercrossingGenerations)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			std::vector<double>::iterator longestPath = std::max_element(pathLengths1.begin(), pathLengths1.end());
			int longestIndex = (int)std::distance(pathLengths1.begin(), longestPath);
			for(int i = 0; i < end - start; i++)
			{
				results(finalCounter, i+start) = intermediate1(longestIndex, i);
			}
		}
	}
	void applyFunnel(int start, int end, int finalCounter, int funnelID, int selfingGenerations)
	{
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHetrozygoteMissingProb = log(hetrozygoteMissingProb);
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerValue = recodedFinals(finalCounter, start);
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
		
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
			{
				int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
				int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2]);
				intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
				if(markerEncodingTheseFounders == startMarkerValue || (startMarkerValue == NA_INTEGER && ((recodedFounders(funnel[founderCounter2], start) == recodedFounders(funnel[founderCounter], start) && homozygoteMissingProb != 0) || (recodedFounders(funnel[founderCounter2], start) != recodedFounders(funnel[founderCounter], start) && hetrozygoteMissingProb != 0))))
				{
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				}
				else
				{
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
			}
		}
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			int markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& previousMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter]];
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter + 1]];
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingMarker = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2]);
					if(encodingMarker == markerValue || (markerValue == NA_INTEGER && ((recodedFounders(funnel[founderCounter2], markerCounter) == recodedFounders(funnel[founderCounter], markerCounter) && homozygoteMissingProb != 0) || (recodedFounders(funnel[founderCounter2], markerCounter) != recodedFounders(funnel[founderCounter], markerCounter) && hetrozygoteMissingProb != 0))))
					{
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousMarker = previousMarkerData.hetData(funnel[founderPreviousCounter], funnel[founderPreviousCounter2]);
								int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2]);
								if(encodingPreviousMarker == previousMarkerValue || (previousMarkerValue == NA_INTEGER && ((recodedFounders(funnel[founderPreviousCounter2], markerCounter) == recodedFounders(funnel[founderPreviousCounter], markerCounter) && homozygoteMissingProb != 0) || (recodedFounders(funnel[founderPreviousCounter2], markerCounter) != recodedFounders(funnel[founderPreviousCounter], markerCounter) && hetrozygoteMissingProb != 0))))
								{
									double multiple = 0;
									if(founderCounter != founderCounter2) multiple += log(2);
									if(founderPreviousCounter != founderPreviousCounter2) multiple += log(2);
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multiple + funnelHaplotypeProbabilities(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
									if(markerValue == NA_INTEGER)
									{
										if(founderCounter2 == founderCounter)
										{
											working[encodingPreviousTheseFounders] += logHomozygoteMissingProb;
										}
										else
										{
											working[encodingPreviousTheseFounders] += logHetrozygoteMissingProb;
										}
									}
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						//This error is no longer valid, because some states are impossible to ever be in - E.g. hetrozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
						//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(markerCounter - start + 1 - identicalIndex));
						intermediate2(encodingTheseFounders, markerCounter-start+1) = encodingTheseFounders;
						pathLengths2[encodingTheseFounders] = *longest;
					}
					else
					{
						pathLengths2[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
				}
			}
#ifndef NDEBUG
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
#endif
			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(key(funnel[0], funnel[0]), identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2]);
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						intermediate2(founderCounter*nFounders + founderCounter2, identicalIndex) = value;
					}
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		/*int markerValue = recodedFinals(finalCounter, start);
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
							working[founderCounter2] = pathLengths1[founderCounter2] + intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter2][founderCounter];
						}
					}
					//Get the longest one, and check that it's not negative infinity.
					std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
					if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
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
		}*/
	}
};
#endif
