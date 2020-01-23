#ifndef FUNNEL_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
#define FUNNEL_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker<nFounders, maxAlleles, false>
{
private:
	funnelHaplotypeToMarker(){}
public:
	static const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
	typedef typename std::array<double, nDifferentProbs> compressedProbabilitiesType;
	template<bool takeLogs> static void convert(rowMajorMatrix<compressedProbabilitiesType>& haplotypeProbabilities, array2<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		int funnel[nFounders];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & ((std::size_t)15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int table[maxAlleles][maxAlleles][nDifferentProbs];
		memset(table, 0, sizeof(int)*maxAlleles*maxAlleles*nDifferentProbs);
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			//firstFounder is the index of the founder within the current funnel
			for(int firstFounder1 = 0; firstFounder1 < nFounders; firstFounder1++)
			{
				for(int firstFounder2 = 0; firstFounder2 < nFounders; firstFounder2++)
				{
					if(firstMarkerPatternData.hetData(funnel[firstFounder1], funnel[firstFounder2]) == firstMarkerValue)
					{
						int index1 = probabilityData<nFounders>::intermediateAllelesMask[firstFounder1][firstFounder2];
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
							{
								for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
								{
									if(secondMarkerPatternData.hetData(funnel[secondFounder1],funnel[secondFounder2]) == secondMarkerValue)
									{
										int index2 = probabilityData<nFounders>::intermediateAllelesMask[secondFounder1][secondFounder2];
										table[firstMarkerValue][secondMarkerValue][probabilityData<nFounders>::intermediateProbabilitiesMask[index1][index2]]++;
									}
								}
							}
						}
					}
				}
			}
		}
		int nPoints = haplotypeProbabilities.getNRows();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			array2<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, selfingGenerationsIndex);
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
	template<bool takeLogs> static void convert16MarkerAlleles(array2<16>& markerProbabilitiesThisRecomb, compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int funnel[16], const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
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
						int index1 = probabilityData<nFounders>::intermediateAllelesMask[firstFounder1][firstFounder2];
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
							{
								for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
								{
									if(secondMarkerPatternData.hetData(funnel[secondFounder1],funnel[secondFounder2]) == secondMarkerValue)
									{
										int index2 = probabilityData<nFounders>::intermediateAllelesMask[secondFounder1][secondFounder2];
										markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilitiesThisRecomb.values[probabilityData<nFounders>::intermediateProbabilitiesMask[index1][index2]];
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
					if(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] == 0) markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = -std::numeric_limits<double>::infinity();
					else markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
		}
	}

};
#endif
