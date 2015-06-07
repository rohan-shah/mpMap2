#ifndef CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#define CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#include "matrices.hpp"
#include "probabilities.hpp"
#include "intercrossingHaplotypeToMarker.hpp"
#include "funnelHaplotypeToMarker.hpp"
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
	xMajorMatrix<array2<maxAlleles> > perFunnelData;
	xMajorMatrix<array2<maxAlleles> > perAIGenerationData;

	rowMajorMatrix<bool> allowableFunnel;
	rowMajorMatrix<bool> allowableAI;
};
template<int maxAlleles> class allMarkerPairData : protected std::vector<singleMarkerPairData<maxAlleles> >
{
public:
	typedef typename std::vector<singleMarkerPairData<maxAlleles> > parent;
	allMarkerPairData(int nMarkerPatternIDs)
		: parent((std::size_t)((nMarkerPatternIDs * nMarkerPatternIDs - nMarkerPatternIDs)/2 + nMarkerPatternIDs), singleMarkerPairData<maxAlleles> (0, 0, 0, 0))
	{}
	singleMarkerPairData<maxAlleles>& operator()(int markerPattern1ID, int markerPattern2ID)
	{
		if(markerPattern1ID < markerPattern2ID) std::swap(markerPattern1ID, markerPattern2ID);
		//Ensure that nMarkerPattern1ID > nMarkerPattern2ID
		int index = markerPattern1ID * (markerPattern1ID + 1) / 2;
		index += markerPattern2ID;
		return parent::operator[](index);
	}
};
template<int maxAlleles, int nFounders> struct constructLookupTableArgs
{
public:
	constructLookupTableArgs(allMarkerPairData<maxAlleles>& computedContributions, markerPatternsToUniqueValuesArgs& markerPatternData)
		: computedContributions(computedContributions), markerPatternData(markerPatternData)
	{}
	allMarkerPairData<maxAlleles>& computedContributions;
	markerPatternsToUniqueValuesArgs& markerPatternData;
	std::vector<funnelEncoding>* funnelEncodings;
	std::vector<double>* recombinationFractions;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
};
template<int maxAlleles> bool isValid(std::vector<array2<maxAlleles> >& markerProbabilities, int nPoints, int nFirstMarkerAlleles, int nSecondMarkerAlleles)
{
	for(int recombCounter1 = 0; recombCounter1 < nPoints; recombCounter1++)
	{
		for(int recombCounter2 = recombCounter1 + 10; recombCounter2 < nPoints; recombCounter2++)
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
			if(sum < 0.005)
			{
				return false;
			} 
		}
	}
	return true;
}
template<int nFounders, int maxAlleles, bool infiniteSelfing> void constructLookupTable(constructLookupTableArgs<maxAlleles, nFounders>& args)
{
	int nMarkerPatternIDs = (int)args.markerPatternData.allMarkerPatterns.size();
	int nRecombLevels = (int)args.recombinationFractions->size();
	int nDifferentFunnels = (int)args.funnelEncodings->size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations->begin(), args.intercrossingGenerations->end());

	int maxSelfing = *std::max_element(args.selfingGenerations->begin(), args.selfingGenerations->end());
	int minSelfing = *std::min_element(args.selfingGenerations->begin(), args.selfingGenerations->end());
	typedef typename expandedProbabilities<nFounders, infiniteSelfing>::type expandedProbabilitiesType;

	//Only compute the haplotype probabilities once. This is for the no intercrossing case
	rowMajorMatrix<expandedProbabilitiesType> funnelHaplotypeProbabilities(nRecombLevels, maxSelfing-minSelfing+1);
	//In order to determine if a marker combination is informative, we use a much finer numerical grid.
	const int nFinerPoints = 101;
	rowMajorMatrix<expandedProbabilitiesType> finerFunnelHaplotypeProbabilities(nFinerPoints, maxSelfing-minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
		{
			expandedGenotypeProbabilities<nFounders, infiniteSelfing>::noIntercross(funnelHaplotypeProbabilities(recombCounter, selfingGenerations-minSelfing), (*args.recombinationFractions)[recombCounter], selfingGenerations);
		}
		for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
		{
			expandedGenotypeProbabilities<nFounders, infiniteSelfing>::noIntercross(finerFunnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0), selfingGenerations);
		}
	}
	//Similarly for the intercrossing generation haplotype probabilities
	xMajorMatrix<expandedProbabilitiesType> intercrossingHaplotypeProbabilities(nRecombLevels, maxAIGenerations, maxSelfing - minSelfing+1);
	xMajorMatrix<expandedProbabilitiesType> finerIntercrossingHaplotypeProbabilities(nFinerPoints, maxAIGenerations, maxSelfing - minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int aiCounter = 1; aiCounter <= maxAIGenerations; aiCounter++)
		{
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				expandedGenotypeProbabilities<nFounders, infiniteSelfing>::withIntercross(intercrossingHaplotypeProbabilities(recombCounter, aiCounter-1, selfingGenerations), aiCounter, (*args.recombinationFractions)[recombCounter], selfingGenerations);
			}
			for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
			{
				expandedGenotypeProbabilities<nFounders, infiniteSelfing>::withIntercross(finerIntercrossingHaplotypeProbabilities(recombCounter, aiCounter-1, selfingGenerations), aiCounter, 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0), selfingGenerations);
			}
		}
	}
#ifdef USE_OPENMP
	#pragma omp parallel 
#endif
	{
		std::vector<array2<maxAlleles> > markerProbabilities(nFinerPoints);
		//These are inputs, used when we want to convert haplotype probabilities (n * n tables) into marker probabilities (these are much smaller tables typically 2*2)
		funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing> finerFunnelHaplotypeToMarker(finerFunnelHaplotypeProbabilities);
		funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing> normalFunnelHaplotypeToMarker(funnelHaplotypeProbabilities);

		intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing> finerIntercrossingHaplotypeToMarker(finerIntercrossingHaplotypeProbabilities);
		intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing> normalIntercrossingHaplotypeToMarker(intercrossingHaplotypeProbabilities);
		//This next loop is a big chunk of code, but does NOT grow with problem size (number of markers, number of lines). Well, it grows but to some fixed limit, because there are only so many marker patterns. 
#ifdef USE_OPENMP
		#pragma omp for schedule(dynamic, 1)
#endif
		for(int firstPattern = 0; firstPattern < nMarkerPatternIDs; firstPattern++)
		{
			markerData& firstMarkerPatternData = args.markerPatternData.allMarkerPatterns[firstPattern];
			for(int secondPattern = 0; secondPattern < nMarkerPatternIDs; secondPattern++)
			{
				markerData& secondMarkerPatternData = args.markerPatternData.allMarkerPatterns[secondPattern];
				//The data for this pair of markers
				singleMarkerPairData<maxAlleles> thisMarkerPairData(nRecombLevels, nDifferentFunnels, maxAIGenerations, maxSelfing - minSelfing + 1);
				for(int selfingCounter = minSelfing; selfingCounter <= maxSelfing; selfingCounter++)
				{
					//Compute marker probabilities for a finer grid. If me seem to see a repeated probability model (numerically, up to a tolerance), then in that particular situtation this pair of markers is no good
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						finerFunnelHaplotypeToMarker.template convert<false>(&(markerProbabilities[0]), (*args.funnelEncodings)[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter);
						thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues);
					}
					for(int intercrossingGeneration = 1; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
					{
						finerIntercrossingHaplotypeToMarker.template convert<false>(&(markerProbabilities[0]), intercrossingGeneration, firstMarkerPatternData, secondMarkerPatternData, selfingCounter);
						thisMarkerPairData.allowableAI(intercrossingGeneration-1, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues);
					}
					//The next two loops relate to the input recombination fractions
					for(int intercrossingGeneration = 1; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
					{
						array2<maxAlleles>* markerProbabilitiesThisIntercrossing = &(thisMarkerPairData.perAIGenerationData(0, intercrossingGeneration-1, selfingCounter - minSelfing));
						if(thisMarkerPairData.allowableAI(intercrossingGeneration-1, selfingCounter - minSelfing))
						{
							normalIntercrossingHaplotypeToMarker.template convert<true>(markerProbabilitiesThisIntercrossing, intercrossingGeneration, firstMarkerPatternData, secondMarkerPatternData, selfingCounter);
						}
					}
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						array2<maxAlleles>* markerProbabilitiesThisFunnel = &(thisMarkerPairData.perFunnelData(0, funnelCounter, selfingCounter - minSelfing));
						memset(&markerProbabilitiesThisFunnel, 0, sizeof(array2<maxAlleles>));
						if(thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing))
						{
							normalFunnelHaplotypeToMarker.template convert<true>(markerProbabilitiesThisFunnel, (*args.funnelEncodings)[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter);
						}
					}
				}
				args.computedContributions(firstPattern, secondPattern).swap(thisMarkerPairData);
			}
		}
	}
}
#endif