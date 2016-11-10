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
#include "joinMapWithExtra.h"
template<int nFounders> struct viterbiAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
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
	int minAIGenerations;
	int maxAIGenerations;
	int maxAlleles;
	Rcpp::IntegerMatrix key;
	double heterozygoteMissingProb, homozygoteMissingProb, errorProb;
	std::vector<array2<nFounders> >* logIntercrossingSingleLociHaplotypeProbabilities, *intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* logFunnelSingleLociHaplotypeProbabilities, *funnelSingleLociHaplotypeProbabilities;
	const positionData& allPositions;
	viterbiAlgorithm(markerPatternsToUniqueValuesArgs& markerData, int maxChromosomeSize, const positionData& allPositions)
		: intermediate1((nFounders*(nFounders+1))/2, maxChromosomeSize), intermediate2(nFounders*nFounders, maxChromosomeSize), pathLengths1((nFounders*(nFounders+1))/2), pathLengths2((nFounders*(nFounders+1))/2), working((nFounders*(nFounders+1)/2)), logIntercrossingHaplotypeProbabilities(NULL), logFunnelHaplotypeProbabilities(NULL), markerData(markerData), lineFunnelIDs(NULL), lineFunnelEncodings(NULL), intercrossingGenerations(NULL), selfingGenerations(NULL), minSelfingGenerations(-1), maxSelfingGenerations(-1), minAIGenerations(-1), maxAIGenerations(-1), maxAlleles(-1), heterozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), homozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), errorProb(std::numeric_limits<double>::quiet_NaN()), logIntercrossingSingleLociHaplotypeProbabilities(NULL), intercrossingSingleLociHaplotypeProbabilities(NULL), logFunnelSingleLociHaplotypeProbabilities(NULL), funnelSingleLociHaplotypeProbabilities(NULL), allPositions(allPositions)
	{}
	void apply(int startPosition, int endPosition)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || lineFunnelIDs == NULL || lineFunnelEncodings == NULL || intercrossingGenerations == NULL || selfingGenerations == NULL)
		{
			throw std::runtime_error("Internal error");
		}
		if((heterozygoteMissingProb != heterozygoteMissingProb || heterozygoteMissingProb == 0) && (homozygoteMissingProb != homozygoteMissingProb || homozygoteMissingProb == 0))
		{
			throw std::runtime_error("One of heterozygoteMissingProb and homozygoteMissingProb must be non-zero");
		}
		if(errorProb < 0 || errorProb >= 1)
		{
			throw std::runtime_error("Input errorProb must be in [0, 1)");
		}
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
		minAIGenerations = *std::min_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		maxAIGenerations = *std::max_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		minAIGenerations = std::max(minAIGenerations, 1);
		int nFinals = recodedFinals.nrow();

		//If there's not meant to be any missing values, check that first
		if(homozygoteMissingProb == 0 && heterozygoteMissingProb == 0)
		{
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				for(int markerCounter = startPosition; markerCounter < endPosition; markerCounter++)
				{
					if(recodedFinals(finalCounter, markerCounter) == NA_INTEGER)
					{
						throw std::runtime_error("Inputs heterozygoteMissingProb and homozygoteMissingProb imply that missing values are not allowed");
					}
				}
			}
		}
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if((*intercrossingGenerations)[finalCounter] == 0)
			{
				applyFunnel(startPosition, endPosition, finalCounter, (*lineFunnelIDs)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			else
			{
				applyIntercrossing(startPosition, endPosition, finalCounter, (*intercrossingGenerations)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			std::vector<double>::iterator longestPath = std::max_element(pathLengths1.begin(), pathLengths1.end());
			int longestIndex = (int)std::distance(pathLengths1.begin(), longestPath);
			for(int i = 0; i < endPosition - startPosition; i++)
			{
				results(finalCounter, i+startPosition) = intermediate1(longestIndex, i) + 1;
			}
		}
	}
	void applyFunnel(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(errorProb == 0)
		{
			applyFunnelNoError(startPosition, endPosition, finalCounter, funnelID, selfingGenerations);
		}
		else 
		{
			applyFunnelWithError(startPosition, endPosition, finalCounter, funnelID, selfingGenerations);
		}
	}
	void applyFunnelNoError(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(logFunnelHaplotypeProbabilities == NULL || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		//Initialise the algorithm
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());

		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			//If so things are a bit more simple. 
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				}
			}
		}
		else
		{
			int startMarkerValue = recodedFinals(finalCounter, startMarkerIndex);
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					pathLengths1[encodingTheseFounders] = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					if(markerEncodingTheseFounders == startMarkerValue)
					{}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) == recodedFounders(funnel[founderCounter], startMarkerIndex))
					{
						if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) != recodedFounders(funnel[founderCounter], startMarkerIndex))
					{
						if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
				}
			}
		}
		int identicalIndex = 0;
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			//Is the current marker just a placeholder, a location where there isn't any data?
			if(markerIndex == -1)
			{
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						//Account for the fact that each heterozygote is only counted once, so the probabilities are half what they should really be. 
						if(founderCounter != founderCounter2) multipleNextMarker += log(2);

						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									double multiplePreviousMarker = 0;
									if(founderPreviousCounter != founderPreviousCounter2) multiplePreviousMarker += log(2);
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multiplePreviousMarker + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
						//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
						intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
						pathLengths2[encodingTheseFounders] = *longest;
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter,markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingMarker = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						//Account for the fact that each heterozygote is only counted once, so the probabilities are half what they should really be. 
						if(founderCounter != founderCounter2) multipleNextMarker += log(2);
						if(markerValue == encodingMarker)
						{}
						else if(markerValue == NA_INTEGER && founderCounter == founderCounter2)
						{
							if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else if(markerValue == NA_INTEGER && founderCounter != founderCounter2)
						{
							if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
						if(multipleNextMarker != -std::numeric_limits<double>::infinity())
						{
							//Founder at the previous marker. 
							std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
							for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
							{
								for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
								{
									int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
									if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
									{
										double multiplePreviousMarker = 0;
										if(founderPreviousCounter != founderPreviousCounter2) multiplePreviousMarker += log(2);
										working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multiplePreviousMarker + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
									}
								}
							}
							//Get the shortest one, and check that it's not negative infinity.
							std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
							//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
							//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
							int bestPrevious = (int)std::distance(working.begin(), longest);
							
							memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
							intermediate2(encodingTheseFounders,positionCounter-startPosition+1) = encodingTheseFounders;
							pathLengths2[encodingTheseFounders] = *longest;
						}
						else
						{
							pathLengths2[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						}
					}
				}
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them. 
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(positionCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != positionCounter-startPosition + 1)
			{
				int value = intermediate1(key(funnel[0], funnel[0])-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyFunnelWithError(int start, int end, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(funnelSingleLociHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		//Initialise the algorithm
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerValue = recodedFinals(finalCounter, start);
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];

		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		double errorTermStart1 = log((1 - errorProb) + errorProb * 1.0 / (double) startMarkerData.nObservedValues);
		double errorTermStart2 = log(errorProb * 1.0 / (double) startMarkerData.nObservedValues);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			//We don't need to consider the other triangle of the matrix startMarkerData.hetData(founderCounter, founderCounter2)
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
				int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
				intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
				pathLengths1[encodingTheseFounders] = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				if(markerEncodingTheseFounders == startMarkerValue)
				{
					pathLengths1[encodingTheseFounders] += errorTermStart1;
				}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], start) == recodedFounders(funnel[founderCounter], start))
				{
					if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], start) != recodedFounders(funnel[founderCounter], start))
				{
					if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else pathLengths1[encodingTheseFounders] += errorTermStart2;
				pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
			}
		}
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter + 1]];

			//The founder at the next marker
			double errorTermCurrentMarker1 = log((1 - errorProb) + errorProb * 1.0 / (double) currentMarkerData.nObservedValues);
			double errorTermCurrentMarker2 = log(errorProb * 1.0 / (double) currentMarkerData.nObservedValues);
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingMarker = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					double multipleNextMarker = 0;
					if(founderCounter != founderCounter2) multipleNextMarker += log(2);
					if(encodingMarker == markerValue)
					{
						multipleNextMarker += errorTermCurrentMarker1;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerCounter) == recodedFounders(funnel[founderCounter], markerCounter))
					{
						if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerCounter) != recodedFounders(funnel[founderCounter], markerCounter))
					{
						if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else multipleNextMarker += errorTermCurrentMarker2;
					if(multipleNextMarker != -std::numeric_limits<double>::infinity())
					{
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									double multiplePreviousMarker = 0;
									if(founderPreviousCounter != founderPreviousCounter2) multiplePreviousMarker += log(2);
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + multiplePreviousMarker + (*logFunnelHaplotypeProbabilities)(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
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
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them. 
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(key(funnel[0], funnel[0])-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(errorProb == 0)
		{
			applyIntercrossingNoError(start, end, finalCounter, intercrossingGeneration, selfingGenerations);
		}
		else
		{
			applyIntercrossingWithError(start, end, finalCounter, intercrossingGeneration, selfingGenerations);
		}
	}
	void applyIntercrossingNoError(int start, int end, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		//Initialise the algorithm
		int startMarkerValue = recodedFinals(finalCounter, start);
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
		
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int markerEncodingTheseFounders = startMarkerData.hetData(founderCounter, founderCounter2);
				int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
				intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
				pathLengths1[encodingTheseFounders] = (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				if(markerEncodingTheseFounders == startMarkerValue)
				{}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, start) == recodedFounders(founderCounter, start))
				{
					if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, start) != recodedFounders(founderCounter, start))
				{
					if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
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
					int encodingMarker = currentMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					double multipleNextMarker = 0;
					if(founderCounter != founderCounter2) multipleNextMarker += log(2);
					if(encodingMarker == markerValue)
					{}
					else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerCounter) == recodedFounders(founderCounter, markerCounter))
					{
						if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerCounter) != recodedFounders(founderCounter, markerCounter))
					{
						if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else multipleNextMarker = -std::numeric_limits<double>::infinity();
					if(multipleNextMarker != -std::numeric_limits<double>::infinity())
					{
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									double multiplePreviousMarker = 0;
									if(founderPreviousCounter != founderPreviousCounter2) multiplePreviousMarker += log(2);
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multiplePreviousMarker + multipleNextMarker + (*logIntercrossingHaplotypeProbabilities)(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
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
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(key(0,0)-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossingWithError(int start, int end, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(intercrossingSingleLociHaplotypeProbabilities == NULL || logIntercrossingHaplotypeProbabilities == NULL || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		//Initialise the algorithm
		int startMarkerValue = recodedFinals(finalCounter, start);
		::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
		
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());
		double errorTermStart1 = log((1 - errorProb) + errorProb * 1.0 / (double) startMarkerData.nObservedValues);
		double errorTermStart2 = log(errorProb * 1.0 / (double) startMarkerData.nObservedValues);
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
			{
				int markerEncodingTheseFounders = startMarkerData.hetData(founderCounter, founderCounter2);
				int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
				intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
				pathLengths1[encodingTheseFounders] = (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				if(markerEncodingTheseFounders == startMarkerValue)
				{
					pathLengths1[encodingTheseFounders] += errorTermStart1;
				}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, start) == recodedFounders(founderCounter, start))
				{
					if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, start) != recodedFounders(founderCounter, start))
				{
					if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
				}
				else
				{
					pathLengths1[encodingTheseFounders] += errorTermStart2;
				}
				pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
			}
		}
		int identicalIndex = 0;
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& previousMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter]];
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter + 1]];

			double errorTermCurrentMarker1 = log((1 - errorProb) + errorProb * 1.0 / (double) currentMarkerData.nObservedValues);
			double errorTermCurrentMarker2 = log(errorProb * 1.0 / (double) currentMarkerData.nObservedValues);;
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingMarker = currentMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					double multipleNextMarker = 0;
					if(founderCounter != founderCounter2) multipleNextMarker += log(2);
					if(encodingMarker == markerValue)
					{
						multipleNextMarker += errorTermCurrentMarker1;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerCounter) == recodedFounders(founderCounter, markerCounter))
					{
						if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerCounter) != recodedFounders(founderCounter, markerCounter))
					{
						if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
					}
					else multipleNextMarker += errorTermCurrentMarker2;
					if(multipleNextMarker != -std::numeric_limits<double>::infinity())
					{
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									double multiplePreviousMarker = 0;
									if(founderPreviousCounter != founderPreviousCounter2) multiplePreviousMarker += log(2);
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + multiplePreviousMarker + (*logIntercrossingHaplotypeProbabilities)(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
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
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(markerCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != markerCounter-start + 1)
			{
				int value = intermediate1(key(0,0)-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
};
#endif
