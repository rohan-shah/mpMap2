#ifndef VITERBI_FINITE_SELFING_HEADER_GUARD
#define VITERBI_FINITE_SELFING_HEADER_GUARD
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
template<int nFounders> struct viterbiAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	rowMajorMatrix<int> intermediate1, intermediate2;
	rowMajorMatrix<bool> error1, error2;
	Rcpp::IntegerMatrix results;
	Rcpp::IntegerMatrix resultsErrors;
	std::vector<double> pathLengths1, pathLengths2;
	std::vector<double> working;
	xMajorMatrix<expandedProbabilitiesType>* logIntercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>* logFunnelHaplotypeProbabilities;
	xMajorMatrix<expandedProbabilitiesType>* intercrossingHaplotypeProbabilities;
	rowMajorMatrix<expandedProbabilitiesType>* funnelHaplotypeProbabilities;
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
	std::function<void(unsigned long long)> updateProgress;
	viterbiAlgorithm(markerPatternsToUniqueValuesArgs& markerData, int maxChromosomeSize, const positionData& allPositions)
		: intermediate1((nFounders*(nFounders+1))/2, maxChromosomeSize), intermediate2(nFounders*nFounders, maxChromosomeSize), error1((nFounders*(nFounders+1))/2, maxChromosomeSize), error2(nFounders*nFounders, maxChromosomeSize), pathLengths1((nFounders*(nFounders+1))/2), pathLengths2((nFounders*(nFounders+1))/2), working((nFounders*(nFounders+1)/2)), logIntercrossingHaplotypeProbabilities(NULL), logFunnelHaplotypeProbabilities(NULL), markerData(markerData), lineFunnelIDs(NULL), lineFunnelEncodings(NULL), intercrossingGenerations(NULL), selfingGenerations(NULL), minSelfingGenerations(-1), maxSelfingGenerations(-1), minAIGenerations(-1), maxAIGenerations(-1), maxAlleles(-1), heterozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), homozygoteMissingProb(std::numeric_limits<double>::quiet_NaN()), errorProb(std::numeric_limits<double>::quiet_NaN()), logIntercrossingSingleLociHaplotypeProbabilities(NULL), logFunnelSingleLociHaplotypeProbabilities(NULL), allPositions(allPositions)
	{}
	void apply(int startPosition, int endPosition)
	{
		if(logIntercrossingHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || lineFunnelIDs == NULL || lineFunnelEncodings == NULL || intercrossingGenerations == NULL || selfingGenerations == NULL || minAIGenerations == -1 || maxAIGenerations == -1)
		{
			throw std::runtime_error("Internal error");
		}
		if((heterozygoteMissingProb != heterozygoteMissingProb || heterozygoteMissingProb == 0) && (homozygoteMissingProb != homozygoteMissingProb || homozygoteMissingProb == 0))
		{
			throw std::runtime_error("One of heterozygoteMissingProb and homozygoteMissingProb must be non-zero");
		}
		if(errorProb != errorProb || errorProb < 0 || errorProb >= 1)
		{
			throw std::runtime_error("Input errorProb must be in [0, 1)");
		}
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
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
			if(errorProb == 0)
			{
				for(int i = 0; i < endPosition - startPosition; i++)
				{
					resultsErrors(finalCounter, i+startPosition) = false;
				}
			}
			else
			{
				for(int i = 0; i < endPosition - startPosition; i++)
				{
					resultsErrors(finalCounter, i+startPosition) = error1(longestIndex, i);
				}
			}
			updateProgress((unsigned long long)finalCounter);
		}
	}
#include "viterbiFiniteSelfingApplyFunnel.h"
#include "viterbiFiniteSelfingApplyIntercrossing.h"
};
#endif
