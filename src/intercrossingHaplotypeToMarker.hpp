#ifndef INTERCROSSING_HAPLOTYPE_TO_MARKER_HEADER_GUARD
#define INTERCROSSING_HAPLOTYPE_TO_MARKER_HEADER_GUARD
template<int nFounders, int maxAlleles, bool infiniteSelfing> struct intercrossingHaplotypeToMarker;
template<int nFounders, int maxAlleles> struct intercrossingHaplotypeToMarker<nFounders, maxAlleles, true>
{
public:
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
	intercrossingHaplotypeToMarker(xMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	xMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities;
	template<bool takeLogs> void convert(array2<maxAlleles>* markerProbabilities, int intercrossingGeneration, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		int nPoints = haplotypeProbabilities.getSizeX();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			array2<nFounders>& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, intercrossingGeneration-1, selfingGenerationsIndex);
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
				//now take logs of every value in markerProbabilities
				for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
				{
					for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
					{
						markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
					}
				}
			}
		}
	}
};
template<int nFounders, int maxAlleles> struct intercrossingHaplotypeToMarker<nFounders, maxAlleles, false>
{
public:
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	intercrossingHaplotypeToMarker(xMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	xMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities;
	template<bool takeLogs> void convert(array2<maxAlleles>* markerProbabilities, int intercrossingGeneration, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		int nPoints = haplotypeProbabilities.getSizeX();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, intercrossingGeneration-1, selfingGenerationsIndex);
			memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<maxAlleles>));
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int firstFounder1 = 0; firstFounder1 < nFounders; firstFounder1++)
				{
					for(int firstFounder2 = 0; firstFounder2 < nFounders; firstFounder2++)
					{
						if(firstMarkerPatternData.hetData(firstFounder1, firstFounder2) == firstMarkerValue)
						{
							for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
							{
								for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
								{
									for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
									{
										if(secondMarkerPatternData.hetData(secondFounder1, secondFounder2) == secondMarkerValue)
										{
											markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilitiesThisRecomb.values[firstFounder1][firstFounder2][secondFounder1][secondFounder2];
										}
									}
								 }
							}
						}
					}
				}
			}
			if(takeLogs)
			{
				//now take logs of every value in markerProbabilities
				for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
				{
					for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
					{
						markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
					}
				}
			}
		}
	}
};
#endif