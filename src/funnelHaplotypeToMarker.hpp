#ifndef FUNNEL_HAPLOTYPE_TO_MARKER_HEADER_GUARD
#define FUNNEL_HAPLOTYPE_TO_MARKER_HEADER_GUARD
template<int nFounders, int maxAlleles, bool infiniteSelfing> struct funnelHaplotypeToMarker;
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker<nFounders, maxAlleles, true>
{
public:
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
	funnelHaplotypeToMarker(rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	int funnel[16];
	rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities;
	template<bool takeLogs> void convert(array2<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getNRows();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			array2<nFounders>& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, selfingGenerationsIndex);
			memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<maxAlleles>));
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				//firstFounder is the index of the founder within the current funnel
				for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
				{
					if(firstMarkerPatternData.hetData(funnel[firstFounder], funnel[firstFounder]) == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
							{
								 if(secondMarkerPatternData.hetData(funnel[secondFounder], funnel[secondFounder]) == secondMarkerValue)
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
	}
};
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker<nFounders, maxAlleles, false>
{
public:
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	funnelHaplotypeToMarker(rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	int funnel[16];
	rowMajorMatrix<expandedProbabilitiesType>& haplotypeProbabilities;
	template<bool takeLogs> void convert(array2<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getNRows();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			expandedProbabilitiesType& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, selfingGenerationsIndex);
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
						markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
					}
				}
#ifndef NDEBUG
				if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
			}
		}
	}
};
#endif
