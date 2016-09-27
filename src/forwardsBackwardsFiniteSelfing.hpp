#ifndef FORWARDS_BACKWARDS_FINITE_SELFING_HEADER_GUARD
#define FORWARDS_BACKWARDS_FINITE_SELFING_HEADER_GUARD
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
template<int nFounders> struct forwardsBackwardsAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	Rcpp::NumericMatrix results;
	xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities;
	markerPatternsToUniqueValuesArgs& markerData;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minSelfingGenerations;
	int maxSelfingGenerations;
	int minAIGenerations, maxAIGenerations;
	double heterozygoteMissingProb, homozygoteMissingProb;
	Rcpp::IntegerMatrix key;
	rowMajorMatrix<double> forwardProbabilities, backwardProbabilities;
	std::vector<array2<nFounders> >* intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* funnelSingleLociHaplotypeProbabilities;
	forwardsBackwardsAlgorithm(markerPatternsToUniqueValuesArgs& markerData, xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities, rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities, int maxChromosomeSize)
		: intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities), markerData(markerData), forwardProbabilities(nFounders*(nFounders+1)/2, maxChromosomeSize), backwardProbabilities(nFounders*(nFounders+1)/2, maxChromosomeSize)
	{}
	void apply(int start, int end)
	{
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
		minAIGenerations = *std::min_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		maxAIGenerations = *std::max_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		minAIGenerations = std::max(minAIGenerations, 1);
		int nFinals = recodedFinals.nrow();
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
		}
	}
	void applyFunnel(int start, int end, int finalCounter, int funnelID, int selfingGenerations)
	{
		//Compute forward probabilities
		int markerValue = recodedFinals(finalCounter, start);
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		{
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					if(markerValue == markerEncodingTheseFounders)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], start) == recodedFounders(funnel[founderCounter], start) && homozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * homozygoteMissingProb;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], start) != recodedFounders(funnel[founderCounter], start) && heterozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * heterozygoteMissingProb;
					}
					else forwardProbabilities(encodingTheseFounders, 0) = 0;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				forwardProbabilities(counter, 0) /= sum;
			}
		}
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter+1]];
			double sum = 0;
			//The founders at the new marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					forwardProbabilities(encodingTheseFounders, markerCounter - start + 1) = 0;
					bool markerMatches = markerValue == markerEncodingTheseFounders;
					bool missingHet = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerCounter + 1) != recodedFounders(funnel[founderCounter], markerCounter + 1) && heterozygoteMissingProb != 0;
					bool missingHomo = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerCounter + 1) == recodedFounders(funnel[founderCounter], markerCounter + 1) && homozygoteMissingProb != 0;
					if(markerMatches || missingHet || missingHomo)
					{
						double factor = 1;
						if(missingHet) factor = heterozygoteMissingProb;
						if(missingHomo) factor = homozygoteMissingProb;
						//Founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								forwardProbabilities(encodingTheseFounders, markerCounter - start + 1) += forwardProbabilities(encodingPreviousFounders, markerCounter - start) * funnelHaplotypeProbabilities(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * factor;
							}
						}
					}
					sum += forwardProbabilities(encodingTheseFounders, markerCounter - start + 1);
				}
			}
			for(int counter = 0; counter < nFounders*(nFounders+1)/2; counter++)
			{
				forwardProbabilities(counter, markerCounter - start + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
				backwardProbabilities(encodingTheseFounders, end - start - 1) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
			}
		}
		for(int markerCounter = end - 2; markerCounter >= start; markerCounter--)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter+1]];
			double sum = 0;
			//The founder at the current marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					backwardProbabilities(encodingTheseFounders, markerCounter - start) = 0;
					//The founders at the previous marker
					for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
					{
						for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
						{
							int markerEncodingPreviousFounders = currentMarkerData.hetData(funnel[founderCounterPrevious], funnel[founderCounterPrevious2]);
							int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
							if(markerValue == markerEncodingPreviousFounders)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * funnelHaplotypeProbabilities(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
							else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerCounter + 1) == recodedFounders(founderCounterPrevious, markerCounter + 1) && homozygoteMissingProb != 0)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * funnelHaplotypeProbabilities(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * homozygoteMissingProb;
							}
							else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerCounter + 1) != recodedFounders(founderCounterPrevious, markerCounter + 1) && heterozygoteMissingProb != 0)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * funnelHaplotypeProbabilities(markerCounter-start, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * heterozygoteMissingProb;
							}
						}
					}
					sum += backwardProbabilities(encodingTheseFounders, markerCounter - start);
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				backwardProbabilities(counter, markerCounter - start) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int markerCounter = start; markerCounter < end; markerCounter++)
		{
			double sum = 0;
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, markerCounter) = backwardProbabilities(counter, markerCounter - start) * forwardProbabilities(counter, markerCounter - start);
				sum += results(((nFounders*(nFounders+1))/2)*finalCounter + counter, markerCounter);
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, markerCounter) /= sum;
			}
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		//Compute forward probabilities
		int markerValue = recodedFinals(finalCounter, start);
		{
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[start]];
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					if(markerValue == markerEncodingTheseFounders)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*intercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					}
					else if(markerValue == NA_INTEGER && recodedFounders((std::size_t)founderCounter2, (std::size_t)start) == recodedFounders((std::size_t)founderCounter, (std::size_t)start) && homozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*intercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * homozygoteMissingProb;
					}
					else if(markerValue == NA_INTEGER && recodedFounders((std::size_t)founderCounter2, (std::size_t)start) != recodedFounders((std::size_t)founderCounter, (std::size_t)start) && heterozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*intercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * heterozygoteMissingProb;
					}
					else forwardProbabilities(encodingTheseFounders, 0) = 0;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				forwardProbabilities(counter, 0) /= sum;
			}
		}
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter+1]];
			double sum = 0;
			//The founders at the new marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = currentMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					forwardProbabilities(encodingTheseFounders, markerCounter - start + 1) = 0;
					bool markerMatches = markerValue == markerEncodingTheseFounders;
					//These rather stupid std::size_t casts are necessary to prevent the comma as being interpreted as a comma, and then calling operator()(int), which returns a whole column. 
					bool missingHet = (markerValue == NA_INTEGER) && (recodedFounders((std::size_t)founderCounter2, (std::size_t)(markerCounter + 1)) != recodedFounders((std::size_t)founderCounter, (std::size_t)(markerCounter + 1))) && (heterozygoteMissingProb != 0);
					bool missingHomo = (markerValue == NA_INTEGER) && (recodedFounders((std::size_t)founderCounter2, (std::size_t)(markerCounter + 1)) == recodedFounders((std::size_t)founderCounter, (std::size_t)(markerCounter + 1))) && (homozygoteMissingProb != 0);
					if(markerMatches || missingHet || missingHomo)
					{
						double factor = 1;
						if(missingHet) factor = heterozygoteMissingProb;
						if(missingHomo) factor = homozygoteMissingProb;
						//Founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(founderCounterPrevious, founderCounterPrevious2)-1;
								forwardProbabilities(encodingTheseFounders, markerCounter - start + 1) += forwardProbabilities(encodingPreviousFounders, markerCounter - start) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * factor;
							}
						}
					}
					sum += forwardProbabilities(encodingTheseFounders, markerCounter - start + 1);
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				forwardProbabilities(counter, markerCounter - start + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
				backwardProbabilities(encodingTheseFounders, end - start - 1) = (*intercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
			}
		}
		for(int markerCounter = end - 2; markerCounter >= start; markerCounter--)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerCounter+1]];
			double sum = 0;
			//The founder at the current marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = currentMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					backwardProbabilities(encodingTheseFounders, markerCounter - start) = 0;
					//The founders at the previous marker
					for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
					{
						for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
						{
							int markerEncodingPreviousFounders = currentMarkerData.hetData(founderCounterPrevious, founderCounterPrevious2);
							int encodingPreviousFounders = key(founderCounterPrevious, founderCounterPrevious2)-1;
							if(markerValue == markerEncodingPreviousFounders)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
							else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerCounter + 1) == recodedFounders(founderCounterPrevious, markerCounter + 1) && homozygoteMissingProb != 0)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * homozygoteMissingProb;
							}
							else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerCounter + 1) != recodedFounders(founderCounterPrevious, markerCounter + 1) && heterozygoteMissingProb != 0)
							{
								backwardProbabilities(encodingTheseFounders, markerCounter - start) += backwardProbabilities(encodingPreviousFounders, markerCounter - start + 1) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * heterozygoteMissingProb;
							}
						}
					}
					sum += backwardProbabilities(encodingTheseFounders, markerCounter - start);
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders + 1))/2; counter++)
			{
				backwardProbabilities(counter, markerCounter - start) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int markerCounter = start; markerCounter < end; markerCounter++)
		{
			double sum = 0;
			for(int counter = 0; counter < (nFounders*(nFounders + 1))/2; counter++)
			{
				results((nFounders*(nFounders+1)/2)*finalCounter + counter, markerCounter) = backwardProbabilities(counter, markerCounter - start) * forwardProbabilities(counter, markerCounter - start);
				sum += results(((nFounders*(nFounders+1))/2)*finalCounter + counter, markerCounter);
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, markerCounter) /= sum;
			}
		}
	}
};
#endif
