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
#include "joinMapWithExtra.h"
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
	int minAIGenerations, maxAIGenerations;
	double heterozygoteMissingProb, homozygoteMissingProb;
	Rcpp::IntegerMatrix key;
	rowMajorMatrix<double> forwardProbabilities, backwardProbabilities;
	std::vector<array2<nFounders> >* intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* funnelSingleLociHaplotypeProbabilities;
	const positionData& allPositions;
	forwardsBackwardsAlgorithm(markerPatternsToUniqueValuesArgs& markerData, xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities, rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities, int maxChromosomeSize, const positionData& allPositions)
		: intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities), markerData(markerData), forwardProbabilities(nFounders, maxChromosomeSize), backwardProbabilities(nFounders, maxChromosomeSize), allPositions(allPositions)
	{}
	void apply(int startPosition, int endPosition)
	{
		minAIGenerations = *std::min_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		maxAIGenerations = *std::max_element(intercrossingGenerations->begin(), intercrossingGenerations->end());
		minAIGenerations = std::max(minAIGenerations, 1);
		int nFinals = recodedFinals.nrow();
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if((*intercrossingGenerations)[finalCounter] == 0)
			{
				applyFunnel(startPosition, endPosition, finalCounter, (*lineFunnelIDs)[finalCounter]);
			}
			else
			{
				applyIntercrossing(startPosition, endPosition, finalCounter, (*intercrossingGenerations)[finalCounter]);
			}
		}
	}
#include "forwardsBackwardsInfiniteSelfingApplyFunnel.hpp"
#include "forwardsBackwardsInfiniteSelfingApplyIntercrossing.hpp"
};
#endif
