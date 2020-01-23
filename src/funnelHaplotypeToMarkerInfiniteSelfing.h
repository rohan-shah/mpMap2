#ifndef FUNNEL_HAPLOTYPE_TO_MARKER_INFINITE_SELFING_HEADER_GUARD
#define FUNNEL_HAPLOTYPE_TO_MARKER_INFINITE_SELFING_HEADER_GUARD
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker<nFounders, maxAlleles, true>
{
private:
	funnelHaplotypeToMarker()
	{}
public:
	static const int nDifferentProbs = compressedProbabilities<nFounders, true>::nDifferentProbs;
	typedef typename std::array<double, nDifferentProbs> compressedProbabilitiesType;
	template<bool takeLogs> static void convert(rowMajorMatrix<compressedProbabilitiesType>& haplotypeProbabilities, array2<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerationsIndex)
	{
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & ((std::size_t)15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getNRows();
		int table[maxAlleles][maxAlleles][nDifferentProbs];
		memset(table, 0, sizeof(int)*maxAlleles*maxAlleles*nDifferentProbs);
		for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
		{
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
					if(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] == 0) markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = -std::numeric_limits<double>::infinity();
					else markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisRecomb.values[firstMarkerValue][secondMarkerValue]);
				}
			}
#ifdef INTERNAL_CHECKS
			if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Joint marker probabilities didn't sum to 1");
#endif
		}
	}

};
#endif
