#ifndef INTERCROSSING_HAPLOTYPE_TO_MARKER_INFINITE_SELFING_HEADER_GUARD
#define INTERCROSSING_HAPLOTYPE_TO_MARKER_INFINITE_SELFING_HEADER_GUARD
#include "array2.h"
template<int nFounders, int maxAlleles> struct intercrossingHaplotypeToMarker<nFounders, maxAlleles, true>
{
private:
	intercrossingHaplotypeToMarker()
	{}
public:
	static const int nDifferentProbs = compressedProbabilities<nFounders, true>::nDifferentProbs;
	typedef typename std::array<double, nDifferentProbs> compressedProbabilitiesType;
	template<bool takeLogs> static void convert(xMajorMatrix<compressedProbabilitiesType>& haplotypeProbabilities, array2<maxAlleles>* markerProbabilities, int intercrossingGenerationsIndex, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex, funnelEncoding enc)
	{
		int nPoints = haplotypeProbabilities.getSizeX();
		int table[maxAlleles][maxAlleles][nDifferentProbs];
		memset(table, 0, sizeof(int)*maxAlleles*maxAlleles*nDifferentProbs);
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
			{
				if(firstMarkerPatternData.hetData(firstFounder, firstFounder) == firstMarkerValue)
				{
					for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
					{
						for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
						{
							if(secondMarkerPatternData.hetData(secondFounder, secondFounder) == secondMarkerValue)
							{
								table[firstMarkerValue][secondMarkerValue][probabilityData<nFounders>::infiniteMask[firstFounder][secondFounder]]++;
							}
						}
					}
				}
			}
		}

		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, intercrossingGenerationsIndex, selfingGenerationsIndex);
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
				{
					double currentMarkerProb = 0;
					int* tableOffset = table[firstMarkerValue][secondMarkerValue];
					for(int differentProbCounter = 0; differentProbCounter < nDifferentProbs; differentProbCounter++)
					{
						currentMarkerProb += tableOffset[differentProbCounter]*haplotypeProbabilitiesThisRecomb[differentProbCounter];
					}
					if(takeLogs)
					{
						if(currentMarkerProb == 0) currentMarkerProb = -std::numeric_limits<double>::infinity();
						else currentMarkerProb = log10(currentMarkerProb);
					}
					markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = currentMarkerProb;
				}
			}
		}
	}
	template<bool takeLogs> static void convert(array2<maxAlleles>& markerProbabilitiesThisRecomb, compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int intercrossingGenerationsIndex, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex, funnelEncoding enc, int minAIGenerations)
	{
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<maxAlleles>));
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
			{
				if(firstMarkerPatternData.hetData(firstFounder, firstFounder) == firstMarkerValue)
				{
					for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
					{
						for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
						{
							 if(secondMarkerPatternData.hetData(secondFounder, secondFounder) == secondMarkerValue)
							 {
								markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilitiesThisRecomb.values[firstFounder][secondFounder];
							 }
						 }
					}
				}
			}
		}
		if(takeLogs)
		{
#ifdef INTERNAL_CHECKS
			double sum = 0;
#endif
			//now take logs of every value in markerProbabilities
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
				{
#ifdef INTERNAL_CHECKS
					sum += markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue];
#endif
					markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifdef INTERNAL_CHECKS
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}
	template<bool takeLogs> static void convert16MarkerAlleles(array2<16>& markerProbabilitiesThisRecomb, compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int intercrossingGenerationsIndex, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<16>));
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
			{
				if(firstMarkerPatternData.hetData(firstFounder, firstFounder) == firstMarkerValue)
				{
					for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
					{
						for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
						{
							 if(secondMarkerPatternData.hetData(secondFounder, secondFounder) == secondMarkerValue)
							 {
								markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilitiesThisRecomb[probabilityData<nFounders>::infiniteMask[firstFounder][secondFounder]];
							 }
						 }
					}
				}
			}
		}
		if(takeLogs)
		{
#ifdef INTERNAL_CHECKS
			double sum = 0;
#endif
			//now take logs of every value in markerProbabilities
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
				{
#ifdef INTERNAL_CHECKS
					sum += markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue];
#endif
					markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifdef INTERNAL_CHECKS
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}

};
#endif
