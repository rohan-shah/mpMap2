#ifndef INTERCROSSING_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
#define INTERCROSSING_HAPLOTYPE_TO_MARKER_FINITE_SELFING_HEADER_GUARD
template<int nFounders, int maxAlleles> struct intercrossingHaplotypeToMarker<nFounders, maxAlleles, false>
{
private:
	intercrossingHaplotypeToMarker()
	{}
public:
	static const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
	typedef typename std::array<double, nDifferentProbs> compressedProbabilitiesType;
	template<bool takeLogs> static void convert(xMajorMatrix<compressedProbabilitiesType>& haplotypeProbabilities, array2<maxAlleles>* markerProbabilities, int intercrossingGenerationsIndex, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex, funnelEncoding enc)
	{
		/*The funnel input only makes a difference to the result if there is a single funnel - If there was more than one funnel we assumed random funnels, in which case every choice of funnel here gives the same result. */
		int funnel[nFounders];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & ((std::size_t)15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getSizeX();
		int table[maxAlleles][maxAlleles][nDifferentProbs];
		memset(table, 0, sizeof(int)*maxAlleles*maxAlleles*nDifferentProbs);
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
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
									if(secondMarkerPatternData.hetData(funnel[secondFounder1], funnel[secondFounder2]) == secondMarkerValue)
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
	template<bool takeLogs> static void convert16MarkerAlleles(array2<16>& markerProbabilitiesThisRecomb, compressedProbabilitiesType& haplotypeProbabilitiesThisRecomb, int intercrossingGenerationsIndex, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex, funnelEncoding enc)
	{
		int funnel[nFounders];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & ((std::size_t)15 << (4*founderCounter))) >> (4*founderCounter));
		}
		memset(&markerProbabilitiesThisRecomb, 0, sizeof(array2<16>));
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
			for(int firstFounder1 = 0; firstFounder1 < nFounders; firstFounder1++)
			{
				for(int firstFounder2 = 0; firstFounder2 < nFounders; firstFounder2++)
				{
					int index1 = probabilityData<nFounders>::intermediateAllelesMask[firstFounder1][firstFounder2];
					if(firstMarkerPatternData.hetData(funnel[firstFounder1], funnel[firstFounder2]) == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder1 = 0; secondFounder1 < nFounders; secondFounder1++)
							{
								for(int secondFounder2 = 0; secondFounder2 < nFounders; secondFounder2++)
								{
									int index2 = probabilityData<nFounders>::intermediateAllelesMask[secondFounder1][secondFounder2];
									if(secondMarkerPatternData.hetData(funnel[secondFounder1], funnel[secondFounder2]) == secondMarkerValue)
									{
										markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilitiesThisRecomb[probabilityData<nFounders>::intermediateProbabilitiesMask[index1][index2]];
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
