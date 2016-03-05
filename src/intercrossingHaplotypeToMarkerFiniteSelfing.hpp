#ifndef INTERCROSSING_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
#define INTERCROSSING_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
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
			convert<takeLogs>(markerProbabilitiesThisRecomb, haplotypeProbabilitiesThisRecomb, intercrossingGeneration, firstMarkerPatternData, secondMarkerPatternData, selfingGenerationsIndex);
		}
	}
	template<bool takeLogs> static void convert(array2<maxAlleles>& markerProbabilitiesThisRecomb, expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int intercrossingGeneration, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
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
#ifndef NDEBUG
			double sum = 0;
#endif
			//now take logs of every value in markerProbabilities
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
				{
#ifndef NDEBUG
					sum += markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue];
#endif
					markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifndef NDEBUG
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}
	template<bool takeLogs> static void convert16MarkerAlleles(array2<16>& markerProbabilitiesThisRecomb, expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int intercrossingGeneration, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<16>));
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
#ifndef NDEBUG
			double sum = 0;
#endif
			//now take logs of every value in markerProbabilities
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
				{
#ifndef NDEBUG
					sum += markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue];
#endif
					markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifndef NDEBUG
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}

};
#endif
