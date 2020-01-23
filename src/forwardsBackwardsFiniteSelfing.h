#ifndef FORWARDS_BACKWARDS_FINITE_SELFING_HEADER_GUARD
#define FORWARDS_BACKWARDS_FINITE_SELFING_HEADER_GUARD
#include "intercrossingAndSelfingGenerations.h"
#include "recodeFoundersFinalsHets.h"
#include "matrices.h"
#include "probabilities.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#include "funnelsToUniqueValues.h"
#include "estimateRFCheckFunnels.h"
#include "markerPatternsToUniqueValues.h"
#include "intercrossingHaplotypeToMarker.h"
#include "funnelHaplotypeToMarker.h"
#include <limits>
#include "joinMapWithExtra.h"
#include "matrix.h"
template<int nFounders> struct forwardsBackwardsAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	integerMatrix recodedFounders, recodedFinals;
	numericMatrix results;
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
	double heterozygoteMissingProb, homozygoteMissingProb, errorProb;
	integerMatrix key;
	rowMajorMatrix<double> forwardProbabilities, backwardProbabilities;
	std::vector<array2<nFounders> >* intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* funnelSingleLociHaplotypeProbabilities;
	const positionData& allPositions;
	forwardsBackwardsAlgorithm(markerPatternsToUniqueValuesArgs& markerData, xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities, rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities, int maxChromosomeSize, const positionData& allPositions)
		: intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities), markerData(markerData), lineFunnelIDs(NULL), intercrossingGenerations(NULL), selfingGenerations(NULL), minSelfingGenerations(-1), maxSelfingGenerations(-1), minAIGenerations(-1), maxAIGenerations(-1), heterozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), homozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), errorProb(std::numeric_limits<double>::quiet_NaN()), forwardProbabilities(nFounders*(nFounders+1)/2, maxChromosomeSize), backwardProbabilities(nFounders*(nFounders+1)/2, maxChromosomeSize), intercrossingSingleLociHaplotypeProbabilities(NULL), funnelSingleLociHaplotypeProbabilities(NULL), allPositions(allPositions)
	{}
	void apply(int startPosition, int endPosition)
	{
		if(errorProb != errorProb || errorProb < 0 || errorProb >= 1)
		{
			throw std::runtime_error("Input errorProb must be in [0, 1)");
		}
		if(minAIGenerations == -1 || maxAIGenerations == -1)
		{
			throw std::runtime_error("Internal error");
		}
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
		int nFinals = recodedFinals.nRow;
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
		}
	}
#include "forwardsBackwardsFiniteSelfingApplyFunnel.h"
#include "forwardsBackwardsFiniteSelfingApplyIntercrossing.h"
};
#endif
