#ifndef CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#define CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#define N_FINER_POINTS 101
#include "matrices.h"
#include "probabilities.h"
#include "intercrossingHaplotypeToMarker.h"
#include "funnelHaplotypeToMarker.h"
#include "getMinAIGenerations.h"
template<int maxAlleles> struct singleMarkerPairData
{
public:
	singleMarkerPairData(int nRecombLevels, int nDifferentFunnels, int nDifferentAIGenerations, int nDifferentSelfingGenerations)
		: perFunnelData(nRecombLevels, nDifferentFunnels, nDifferentSelfingGenerations), perAIGenerationData(nRecombLevels, nDifferentAIGenerations, nDifferentSelfingGenerations), allowableFunnel(nDifferentFunnels, nDifferentSelfingGenerations), allowableAI(nDifferentAIGenerations, nDifferentSelfingGenerations)
	{}
	void swap(singleMarkerPairData<maxAlleles>& other)
	{
		perFunnelData.swap(other.perFunnelData);
		perAIGenerationData.swap(other.perAIGenerationData);
		allowableFunnel.swap(other.allowableFunnel);
		allowableAI.swap(other.allowableAI);
	}
	singleMarkerPairData()
	{}
	singleMarkerPairData(singleMarkerPairData<maxAlleles>&& other)
		:perFunnelData(std::move(other.perFunnelData)), perAIGenerationData(std::move(other.perAIGenerationData)), allowableFunnel(std::move(other.allowableFunnel)), allowableAI(std::move(other.allowableAI))
	{}
	singleMarkerPairData<maxAlleles>& operator=(singleMarkerPairData<maxAlleles>&& other)
	{
		perFunnelData = std::move(other.perFunnelData);
		perAIGenerationData = std::move(other.perAIGenerationData);
		allowableFunnel = std::move(other.allowableFunnel);
		allowableAI = std::move(other.allowableAI);
	}
	xMajorMatrix<array2<maxAlleles> > perFunnelData;
	xMajorMatrix<array2<maxAlleles> > perAIGenerationData;

	rowMajorMatrix<bool> allowableFunnel;
	rowMajorMatrix<bool> allowableAI;
private:
	singleMarkerPairData(const singleMarkerPairData<maxAlleles>& other);
	singleMarkerPairData<maxAlleles>& operator=(const singleMarkerPairData<maxAlleles>& other);
};
template<int maxAlleles> class allMarkerPairData : protected std::vector<singleMarkerPairData<maxAlleles>* >
{
public:
	typedef typename std::vector<singleMarkerPairData<maxAlleles>* > parent;
	allMarkerPairData(int nMarkerPatternIDs)
		: parent((std::size_t)((nMarkerPatternIDs * nMarkerPatternIDs - nMarkerPatternIDs)/2 + nMarkerPatternIDs), NULL)
	{}
	~allMarkerPairData()
	{
		for(typename parent::iterator i = parent::begin(); i != parent::end(); i++)
		{
			if(*i != NULL)
			{
				delete *i;
			}
		}
		parent::clear();
	}
	singleMarkerPairData<maxAlleles>*& operator()(int markerPattern1ID, int markerPattern2ID)
	{
		if(markerPattern1ID < markerPattern2ID) std::swap(markerPattern1ID, markerPattern2ID);
		//Ensure that nMarkerPattern1ID > nMarkerPattern2ID
		int index = markerPattern1ID * (markerPattern1ID + 1) / 2;
		index += markerPattern2ID;
		return parent::operator[](index);
	}
private:
	allMarkerPairData();
	allMarkerPairData(const allMarkerPairData& other);
	allMarkerPairData& operator=(const allMarkerPairData& other);
};
template<int maxAlleles> struct constructLookupTableArgs
{
public:
	constructLookupTableArgs(allMarkerPairData<maxAlleles>& computedContributions, markerPatternsToUniqueValuesArgs& markerPatternData, const std::vector<int>& rowPatterns, const std::vector<int>& columnPatterns)
		: computedContributions(computedContributions), markerPatternData(markerPatternData), rowPatterns(rowPatterns), columnPatterns(columnPatterns)
	{}
	allMarkerPairData<maxAlleles>& computedContributions;
	markerPatternsToUniqueValuesArgs& markerPatternData;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<funnelEncoding>* allFunnelEncodings;

	const std::vector<double>* recombinationFractions;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	const std::vector<int>& rowPatterns;
	const std::vector<int>& columnPatterns;
};
template<int maxAlleles> bool isValid(std::vector<array2<maxAlleles> >& markerProbabilities, int nPoints, int nFirstMarkerAlleles, int nSecondMarkerAlleles, std::vector<double>& recombLevels)
{
	for(int recombCounter1 = 0; recombCounter1 < nPoints; recombCounter1++)
	{
		for(int recombCounter2 = recombCounter1; recombCounter2 < nPoints; recombCounter2++)
		{
			if(fabs(recombLevels[recombCounter1] - recombLevels[recombCounter2]) > 0.06)
			{
				double sum = 0;
				for(int i = 0; i < nFirstMarkerAlleles; i++)
				{
					for(int j = 0; j < nSecondMarkerAlleles; j++)
					{
						sum += fabs(markerProbabilities[recombCounter1].values[i][j] - markerProbabilities[recombCounter2].values[i][j]);
					}
				}
				//If two different recombination fractions give similar models then this pair of markers is not good
				if(sum < 0.003)
				{
					return false;
				}
			}
		}
	}
	return true;
}
template<int nFounders, int maxAlleles, bool infiniteSelfing> void constructLookupTable(constructLookupTableArgs<maxAlleles>& args)
{
	int nRecombLevels = (int)args.recombinationFractions->size();
	int nDifferentFunnels = (int)args.lineFunnelEncodings->size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations->begin(), args.intercrossingGenerations->end());
	int minAIGenerations = getMinAIGenerations(args.intercrossingGenerations);

	int maxSelfing = *std::max_element(args.selfingGenerations->begin(), args.selfingGenerations->end());
	int minSelfing = *std::min_element(args.selfingGenerations->begin(), args.selfingGenerations->end());

	//Only compute the compressed haplotype probabilities once. This is for the no intercrossing case
	typedef std::array<double, compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs> compressedProbabilitiesType;
	rowMajorMatrix<compressedProbabilitiesType> funnelHaplotypeProbabilities(nRecombLevels, maxSelfing - minSelfing + 1);
	//In order to determine if a marker combination is informative, we use a much finer numerical grid.
	const int nFinerPoints = N_FINER_POINTS;
	std::vector<double> finerRecombLevels(nFinerPoints);
	for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
	{
		finerRecombLevels[recombCounter] = 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0);
	}

	rowMajorMatrix<compressedProbabilitiesType> finerFunnelHaplotypeProbabilities(nFinerPoints, maxSelfing-minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(funnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), (*args.recombinationFractions)[recombCounter], selfingGenerations, args.allFunnelEncodings->size());
		}
		for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(finerFunnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), finerRecombLevels[recombCounter], selfingGenerations, args.allFunnelEncodings->size());
		}
	}
	//Similarly for the intercrossing generation haplotype probabilities
	xMajorMatrix<compressedProbabilitiesType> intercrossingHaplotypeProbabilities(nRecombLevels, maxAIGenerations - minAIGenerations + 1, maxSelfing - minSelfing+1);
	xMajorMatrix<compressedProbabilitiesType> finerIntercrossingHaplotypeProbabilities(nFinerPoints, maxAIGenerations - minAIGenerations + 1, maxSelfing - minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int aiCounter = minAIGenerations; aiCounter <= maxAIGenerations; aiCounter++)
		{
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(intercrossingHaplotypeProbabilities(recombCounter, aiCounter-minAIGenerations, selfingGenerations-minSelfing), aiCounter, (*args.recombinationFractions)[recombCounter], selfingGenerations, args.allFunnelEncodings->size());
			}
			for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(finerIntercrossingHaplotypeProbabilities(recombCounter, aiCounter - minAIGenerations, selfingGenerations - minSelfing), aiCounter, finerRecombLevels[recombCounter], selfingGenerations, args.allFunnelEncodings->size());
			}
		}
	}
#ifdef USE_OPENMP
	#pragma omp parallel 
#endif
	{
		std::vector<array2<maxAlleles> > markerProbabilities(nFinerPoints);
		//This next loop is a big chunk of code, but does NOT grow with problem size (number of markers, number of lines). Well, it grows but to some fixed limit, because there are only so many marker patterns. 
#ifdef USE_OPENMP
#pragma omp for schedule(dynamic)
#endif
		for(std::vector<int>::const_iterator rowPatternIterator = args.rowPatterns.begin(); rowPatternIterator < args.rowPatterns.end(); rowPatternIterator++)
		{
			int firstPattern = *rowPatternIterator;
			for(std::vector<int>::const_iterator columnPatternIterator = args.columnPatterns.begin(); columnPatternIterator != args.columnPatterns.end(); columnPatternIterator++)
			{
				int secondPattern = *columnPatternIterator;
				int copiedFirstPattern = firstPattern;
				int copiedSecondPattern = secondPattern;
				if(copiedSecondPattern < copiedFirstPattern)
				{
					std::swap(copiedFirstPattern, copiedSecondPattern);
				}
				markerData& firstMarkerPatternData = args.markerPatternData.allMarkerPatterns[copiedFirstPattern];
				markerData& secondMarkerPatternData = args.markerPatternData.allMarkerPatterns[copiedSecondPattern];
				//The data for this pair of markers
				if(args.computedContributions(copiedFirstPattern, copiedSecondPattern) != NULL) continue;
				args.computedContributions(copiedFirstPattern, copiedSecondPattern) = new singleMarkerPairData<maxAlleles>(nRecombLevels, nDifferentFunnels, std::max(maxAIGenerations-minAIGenerations+1, 0), maxSelfing - minSelfing + 1);
				singleMarkerPairData<maxAlleles>& thisMarkerPairData = *(args.computedContributions(copiedFirstPattern, copiedSecondPattern));
				for(int selfingCounter = minSelfing; selfingCounter <= maxSelfing; selfingCounter++)
				{
					//Compute marker probabilities for a finer grid. If me seem to see a repeated probability model (numerically, up to a tolerance), then in that particular situtation this pair of markers is no good
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<false>(finerFunnelHaplotypeProbabilities, &(markerProbabilities[0]), (*args.lineFunnelEncodings)[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing);
						thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues, finerRecombLevels);
					}
					for(int intercrossingGeneration = minAIGenerations; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
					{
						intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<false>(finerIntercrossingHaplotypeProbabilities, &(markerProbabilities[0]), intercrossingGeneration - minAIGenerations, firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing, (*args.allFunnelEncodings)[0]);
						thisMarkerPairData.allowableAI(intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues, finerRecombLevels);
					}
					//The next two loops relate to the input recombination fractions
					for(int intercrossingGeneration = minAIGenerations; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
					{
						array2<maxAlleles>* markerProbabilitiesThisIntercrossing = &(thisMarkerPairData.perAIGenerationData(0, intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing));
						if(thisMarkerPairData.allowableAI(intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing))
						{
							intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<true>(intercrossingHaplotypeProbabilities, markerProbabilitiesThisIntercrossing, intercrossingGeneration-minAIGenerations, firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing, (*args.allFunnelEncodings)[0]);
						}
					}
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						array2<maxAlleles>* markerProbabilitiesThisFunnel = &(thisMarkerPairData.perFunnelData(0, funnelCounter, selfingCounter - minSelfing));
						memset(markerProbabilitiesThisFunnel->values, 0, sizeof(array2<maxAlleles>::values));
						if(thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing))
						{
							funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<true>(funnelHaplotypeProbabilities, markerProbabilitiesThisFunnel, (*args.lineFunnelEncodings)[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing);
						}
					}
				}
			}
		}
	}
}
#endif
