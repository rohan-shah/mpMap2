#ifndef FORWARDS_BACKWARDS_INFINITE_SELFING_HEADER_GUARD
#define FORWARDS_BACKWARDS_INFINITE_SELFING_HEADER_GUARD
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
template<int nFounders> struct forwardsBackwardsAlgorithm<nFounders, true>
{
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
	integerMatrix recodedFounders, recodedFinals;
	numericMatrix results;
	xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities;
	markerPatternsToUniqueValuesArgs& markerData;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minAIGenerations, maxAIGenerations;
	double heterozygoteMissingProb, homozygoteMissingProb, errorProb;
	integerMatrix key;
	rowMajorMatrix<double> forwardProbabilities, backwardProbabilities;
	std::vector<array2<nFounders> >* intercrossingSingleLociHaplotypeProbabilities;
	std::vector<array2<nFounders> >* funnelSingleLociHaplotypeProbabilities;
	const positionData& allPositions;
	forwardsBackwardsAlgorithm(markerPatternsToUniqueValuesArgs& markerData, xMajorMatrix<expandedProbabilitiesType>& intercrossingHaplotypeProbabilities, rowMajorMatrix<expandedProbabilitiesType>& funnelHaplotypeProbabilities, int maxChromosomeSize, const positionData& allPositions)
		: intercrossingHaplotypeProbabilities(intercrossingHaplotypeProbabilities), funnelHaplotypeProbabilities(funnelHaplotypeProbabilities), markerData(markerData), lineFunnelIDs(NULL), lineFunnelEncodings(NULL), intercrossingGenerations(NULL), selfingGenerations(NULL), minAIGenerations(-1), maxAIGenerations(-1), heterozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), homozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), errorProb(std::numeric_limits<double>::quiet_NaN()), forwardProbabilities(nFounders, maxChromosomeSize), backwardProbabilities(nFounders, maxChromosomeSize), intercrossingSingleLociHaplotypeProbabilities(NULL), funnelSingleLociHaplotypeProbabilities(NULL), allPositions(allPositions)
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
		int nFinals = recodedFinals.nRow;
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
#include "forwardsBackwardsInfiniteSelfingApplyFunnel.h"
#include "forwardsBackwardsInfiniteSelfingApplyIntercrossing.h"
};
#endif
