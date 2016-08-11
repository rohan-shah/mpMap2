#ifndef FORWARDS_BACKWARDS_INFINITE_SELFING_HEADER_GUARD
#define FORWARDS_BACKWARDS_INFINITE_SELFING_HEADER_GUARD
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
template<int nFounders> struct forwardsBackwardsAlgorithm<nFounders, true>
{
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
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
		: intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities), markerData(markerData), forwardProbabilities(nFounders, maxChromosomeSize), backwardProbabilities(nFounders, maxChromosomeSize)
	{}
	void apply(int start, int end)
	{
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
		}
	}
	void applyFunnel(int start, int end, int finalCounter, int funnelID)
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
			int validInitial = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(funnel[founderCounter], start) == markerValue || markerValue == NA_INTEGER)
				{
					forwardProbabilities(funnel[founderCounter], 0) = 1;
					validInitial++;
				}
				else forwardProbabilities(founderCounter, 0) = 0;
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) /= (double)validInitial;
			}
		}
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			double sum = 0;
			//The founder at the new marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, markerCounter - start + 1) = 0;
				if(recodedFounders(funnel[founderCounter], markerCounter - start + 1) == markerValue || markerValue == NA_INTEGER)
				{
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						forwardProbabilities(funnel[founderCounter], markerCounter - start + 1) += forwardProbabilities(funnel[founderCounter2], markerCounter - start) * funnelHaplotypeProbabilities(markerCounter-start, 0).values[founderCounter2][founderCounter];
					}
				}
				sum += forwardProbabilities(funnel[founderCounter], markerCounter - start + 1);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(funnel[founderCounter], markerCounter - start + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			backwardProbabilities(founderCounter, end - start - 1) = 1/(double)nFounders;
		}
		for(int markerCounter = end - 2; markerCounter >= start; markerCounter--)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			double sum = 0;
			//The founder at the current marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(funnel[founderCounter], markerCounter - start) = 0;
				//The founder at the previous marker
				for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
				{
					if(recodedFounders(funnel[founderCounter2], markerCounter - start + 1) == markerValue || markerValue == NA_INTEGER)
					{
						backwardProbabilities(funnel[founderCounter], markerCounter - start) += backwardProbabilities(funnel[founderCounter2], markerCounter - start + 1) * funnelHaplotypeProbabilities(markerCounter-start, 0).values[founderCounter2][founderCounter];
					}
				}
				sum += backwardProbabilities(funnel[founderCounter], markerCounter - start);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(funnel[founderCounter], markerCounter - start) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int markerCounter = start; markerCounter < end; markerCounter++)
		{
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(4*finalCounter + founderCounter, markerCounter) = backwardProbabilities(founderCounter, markerCounter - start) * forwardProbabilities(founderCounter, markerCounter - start);
				sum += results(4*finalCounter + founderCounter, markerCounter);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(4*finalCounter + founderCounter, markerCounter) /= sum;
			}
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGeneration)
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
			int validInitial = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(funnel[founderCounter], start) == markerValue || markerValue == NA_INTEGER)
				{
					forwardProbabilities(funnel[founderCounter], 0) = 1;
					validInitial++;
				}
				else forwardProbabilities(founderCounter, 0) = 0;
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) /= (double)validInitial;
			}
		}
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			double sum = 0;
			//The founder at the new marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, markerCounter - start + 1) = 0;
				if(recodedFounders(funnel[founderCounter], markerCounter - start + 1) == markerValue || markerValue == NA_INTEGER)
				{
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						forwardProbabilities(funnel[founderCounter], markerCounter - start + 1) += forwardProbabilities(funnel[founderCounter2], markerCounter - start) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
				}
				sum += forwardProbabilities(funnel[founderCounter], markerCounter - start + 1);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(funnel[founderCounter], markerCounter - start + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			backwardProbabilities(founderCounter, end - start - 1) = 1/(double)nFounders;
		}
		for(int markerCounter = end - 2; markerCounter >= start; markerCounter--)
		{
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			double sum = 0;
			//The founder at the current marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(funnel[founderCounter], markerCounter - start) = 0;
				//The founder at the previous marker
				for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
				{
					if(recodedFounders(funnel[founderCounter2], markerCounter - start + 1) == markerValue || markerValue == NA_INTEGER)
					{
						backwardProbabilities(funnel[founderCounter], markerCounter - start) += backwardProbabilities(funnel[founderCounter2], markerCounter - start + 1) * intercrossingHaplotypeProbabilities(markerCounter-start, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
				}
				sum += backwardProbabilities(funnel[founderCounter], markerCounter - start);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(funnel[founderCounter], markerCounter - start) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int markerCounter = start; markerCounter < end; markerCounter++)
		{
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(4*finalCounter + founderCounter, markerCounter) = backwardProbabilities(founderCounter, markerCounter - start) * forwardProbabilities(founderCounter, markerCounter - start);
				sum += results(4*finalCounter + founderCounter, markerCounter);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(4*finalCounter + founderCounter, markerCounter) /= sum;
			}
		}
	}
};
#endif
