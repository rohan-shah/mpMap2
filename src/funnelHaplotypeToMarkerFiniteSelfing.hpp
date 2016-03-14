#ifndef FUNNEL_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
#define FUNNEL_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker<nFounders, maxAlleles, false>
{
public:
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	funnelHaplotypeToMarker(rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities;
	template<bool takeLogs> void convert(array2<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getNRows();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, selfingGenerationsIndex);
			convert<takeLogs>(markerProbabilitiesThisRecomb, haplotypeProbabilitiesThisRecomb, funnel, firstMarkerPatternData, secondMarkerPatternData, selfingGenerationsIndex);
		}
	}
	template<bool takeLogs> static void convert(array2<maxAlleles>& markerProbabilitiesThisRecomb, expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int funnel[16], const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<maxAlleles>));
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			//firstFounder is the index of the founder within the current funnel
			for(int firstFounder1 = 0; firstFounder1 < nFounders; firstFounder1++)
			{
				for(int firstFounder2 = 0; firstFounder2 < nFounders; firstFounder2++)
				{
					if(firstMarkerPatternData.hetData(funnel[firstFounder1], funnel[firstFounder2]) == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
							{
								for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
								{
									if(secondMarkerPatternData.hetData(funnel[secondFounder1],funnel[secondFounder2]) == secondMarkerValue)
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
					if(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] == 0) markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = -std::numeric_limits<double>::infinity();
					else markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifndef NDEBUG
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}
	template<bool takeLogs> static void convert16MarkerAlleles(array2<16>& markerProbabilitiesThisRecomb, expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int funnel[16], const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<16>));
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			//firstFounder is the index of the founder within the current funnel
			for(int firstFounder1 = 0; firstFounder1 < nFounders; firstFounder1++)
			{
				for(int firstFounder2 = 0; firstFounder2 < nFounders; firstFounder2++)
				{
					if(firstMarkerPatternData.hetData(funnel[firstFounder1], funnel[firstFounder2]) == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
							{
								for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
								{
									if(secondMarkerPatternData.hetData(funnel[secondFounder1],funnel[secondFounder2]) == secondMarkerValue)
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
					if(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] == 0) markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = -std::numeric_limits<double>::infinity();
					else markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifndef NDEBUG
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}

};
#endif
