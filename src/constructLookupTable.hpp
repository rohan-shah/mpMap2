#ifndef CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#define CONSTRUCT_LOOKUP_TABLE_HEADER_GUARD
#include "matrices.hpp"
template<int n> struct arrayType
{
public:
	arrayType()
	{}
	double values[n][n];
};
//Note that if we change the mask[8][8] values of 2 to 1 we get mask4 in the first 4x4 block. 
//const int mask4[4][4] = {{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
const int mask[8][8] =
	{
			{0, 1, 2, 2, 2, 2, 2, 2},
			{1, 0, 2, 2, 2, 2, 2, 2},
			{2, 2, 0, 1, 2, 2, 2, 2},
			{2, 2, 1, 0, 2, 2, 2, 2},
			{2, 2, 2, 2, 0, 1, 2, 2},
			{2, 2, 2, 2, 1, 0, 2, 2},
			{2, 2, 2, 2, 2, 2, 0, 1},
			{2, 2, 2, 2, 2, 2, 1, 0}
	};
/*See Karl Bromans paper on intermediate generations. In terms of indices, 
Zero = homozygote, One = other homozygote, Two = hetrozygote
In terms of values, see Table one of the paper. Zero = first equation of table, one = second equation, etc. Note that
we combine equations 4 and 5 into a single state. 
*/
const int rilIntermediateMask[3][3] = 
	{
		{0, 1, 2},
		{1, 0, 2},
		{2, 2, 3}
	};
//Templated function to work out the two-point probabilities with the given recombination fraction (and number of AI generations). Templating allows the number of founders to be a compile-time constant
template<int nFounders> void genotypeProbabilitiesNoIntercross(double (&prob)[3], double recombinationFraction);
template<int nFounders> void genotypeProbabilitiesWithIntercross(double (&prob)[3], int nAIGenarations, double recombinationFraction);
//There are only really two probability values for the 4-way design, but if we put in three values we can use the same mask as for the 8-way case
template<> void genotypeProbabilitiesNoIntercross<4>(double (&prob)[3], double r)
{
	prob[0] = (1-r)/(4+8*r);
	prob[1] = prob[2] = r/(4+8*r);
}
template<> void genotypeProbabilitiesNoIntercross<8>(double (&prob)[3], double r)
{
	prob[0] = (1-r)*(1-r)/(8+16*r);
	prob[1] = r*(1-r)/(8+16*r);
	prob[2] = r/(16+32*r);
}
template<> void genotypeProbabilitiesNoIntercross<2>(double (&prob)[3], double r)
{
	prob[0] = 1/(2*(1 + 2*r));
	prob[1] = r/(1 + 2 * r);
}
template<> void genotypeProbabilitiesWithIntercross<2>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations - 1);
	//calculated by taking the 4-way case and setting both pairs of founders to be identical
	prob[0] = (1/(1 + 2 * r)) * ((1-r)*tmp/2 + (2*r + 1 - tmp) /4);
	prob[1] = (1 - prob[0]*2)/2;
}
template<> void genotypeProbabilitiesWithIntercross<4>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	//prob[0] = (pow(1-r, 1+nAIGenerations)/4+(2*r+1-pow(1-r, nAIGenerations-1))/16)/(1+2*r); 
	prob[0] = (tmp *(1-r)*(1-r)/4 + (2*r + 1 - tmp)/16)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 4 * prob[0]) / 12;
}
template<> void genotypeProbabilitiesWithIntercross<8>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	prob[0] = (tmp *(1-r)*(1-r)*(1-r)/8 + (2*r + 1 - tmp)/64)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 8 * prob[0]) / 56;
}
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesNoIntercross(arrayType<nFounders>& expandedProbabilities, double r, int selfingGenerations)
{
	double probabilities[3];
	genotypeProbabilitiesNoIntercross<nFounders>(probabilities, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities.values[i][j] = probabilities[mask[i][j]];
		}
	}
}
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesWithIntercross(arrayType<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations)
{
	double probabilities[3];
	genotypeProbabilitiesWithIntercross<nFounders>(probabilities, nAIGenerations, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities.values[i][j] = probabilities[mask[i][j]];
		}
	}
}
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
	xMajorMatrix<arrayType<maxAlleles> > perFunnelData;
	xMajorMatrix<arrayType<maxAlleles> > perAIGenerationData;

	rowMajorMatrix<bool> allowableFunnel;
	rowMajorMatrix<bool> allowableAI;
};
template<int maxAlleles> class allMarkerPairData : protected std::vector<singleMarkerPairData<maxAlleles> >
{
public:
	typedef typename std::vector<singleMarkerPairData<maxAlleles> > parent;
	allMarkerPairData(int nMarkerPatternIDs)
		: parent((nMarkerPatternIDs * nMarkerPatternIDs - nMarkerPatternIDs)/2 + nMarkerPatternIDs, singleMarkerPairData<maxAlleles> (0, 0, 0, 0))
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
template<int maxAlleles> struct constructLookupTableArgs
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
template<int maxAlleles> bool isValid(std::vector<arrayType<maxAlleles> >& markerProbabilities, int nPoints, int nFirstMarkerAlleles, int nSecondMarkerAlleles)
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
template<int nFounders, int maxAlleles> struct funnelHaplotypeToMarker
{
public:
	funnelHaplotypeToMarker(rowMajorMatrix<arrayType<nFounders> >& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	int funnel[8];
	rowMajorMatrix<arrayType<nFounders> >& haplotypeProbabilities;
	template<bool takeLogs> void convert(arrayType<maxAlleles>* markerProbabilities, funnelEncoding enc, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerations)
	{
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (7 << (3*founderCounter))) >> (3*founderCounter));
		}
		int nPoints = haplotypeProbabilities.getNRows();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			arrayType<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			arrayType<nFounders>& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, selfingGenerations);
			memset(&markerProbabilitiesThisRecomb, 0, sizeof(arrayType<maxAlleles>));
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				//firstFounder is the index of the founder within the current funnel
				for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
				{
					if(firstMarkerPatternData.founderAlleles[funnel[firstFounder]] == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
							{
								 if(secondMarkerPatternData.founderAlleles[funnel[secondFounder]] == secondMarkerValue)
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
template<int nFounders, int maxAlleles> struct intercrossingHaplotypeToMarker
{
public:
	intercrossingHaplotypeToMarker(xMajorMatrix<arrayType<nFounders> >& haplotypeProbabilities)
		: haplotypeProbabilities(haplotypeProbabilities)
	{}
	xMajorMatrix<arrayType<nFounders> >& haplotypeProbabilities;
	template<bool takeLogs> void convert(arrayType<maxAlleles>* markerProbabilities, int intercrossingGeneration, const markerData& firstMarkerPatternData, const markerData& secondMarkerPatternData, int selfingGenerations)
	{
		int nPoints = haplotypeProbabilities.getSizeX();
		for(int recombCounter = 0; recombCounter < nPoints; recombCounter++)
		{
			arrayType<maxAlleles>& markerProbabilitiesThisRecomb = markerProbabilities[recombCounter];
			arrayType<nFounders>& haplotypeProbabilitiesThisRecomb = haplotypeProbabilities(recombCounter, intercrossingGeneration-1, selfingGenerations);
			memset(&markerProbabilitiesThisRecomb, 0, sizeof(arrayType<maxAlleles>));
			for(int firstMarkerValue = 0; firstMarkerValue < firstMarkerPatternData.nObservedValues; firstMarkerValue++)
			{
				for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
				{
					if(firstMarkerPatternData.founderAlleles[firstFounder] == firstMarkerValue)
					{
						for(int secondMarkerValue = 0; secondMarkerValue < secondMarkerPatternData.nObservedValues; secondMarkerValue++)
						{
							for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
							{
								 if(secondMarkerPatternData.founderAlleles[secondFounder] == secondMarkerValue)
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
template<int nFounders, int maxAlleles, bool infiniteSelfing> void constructLookupTable(constructLookupTableArgs<maxAlleles>& args)
{
	int nMarkerPatternIDs = args.markerPatternData.allMarkerPatterns.size();
	int nRecombLevels = args.recombinationFractions->size();
	int nDifferentFunnels = args.funnelEncodings->size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations->begin(), args.intercrossingGenerations->end());

	int maxSelfing = *std::max_element(args.selfingGenerations->begin(), args.selfingGenerations->end());
	int minSelfing = *std::min_element(args.selfingGenerations->begin(), args.selfingGenerations->end());

	//Only compute the haplotype probabilities once. This is for the no intercrossing case
	rowMajorMatrix<arrayType<nFounders> > funnelHaplotypeProbabilities(nRecombLevels, maxSelfing-minSelfing+1);
	//In order to determine if a marker combination is informative, we use a much finer numerical grid.
	const int nFinerPoints = 101;
	rowMajorMatrix<arrayType<nFounders> > finerFunnelHaplotypeProbabilities(nFinerPoints, maxSelfing-minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(funnelHaplotypeProbabilities(recombCounter, selfingGenerations-minSelfing), (*args.recombinationFractions)[recombCounter], selfingGenerations);
		}
		for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(finerFunnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0), selfingGenerations);
		}
	}
	//Similarly for the intercrossing generation haplotype probabilities
	xMajorMatrix<arrayType<nFounders> > intercrossingHaplotypeProbabilities(nRecombLevels, maxAIGenerations, maxSelfing - minSelfing+1);
	xMajorMatrix<arrayType<nFounders> > finerIntercrossingHaplotypeProbabilities(nFinerPoints, maxAIGenerations, maxSelfing - minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int aiCounter = 1; aiCounter <= maxAIGenerations; aiCounter++)
		{
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(intercrossingHaplotypeProbabilities(recombCounter, aiCounter-1, selfingGenerations), aiCounter, (*args.recombinationFractions)[recombCounter], selfingGenerations);
			}
			for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(finerIntercrossingHaplotypeProbabilities(recombCounter, aiCounter-1, selfingGenerations), aiCounter, 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0), selfingGenerations);
			}
		}
	}
#ifdef USE_OPENMP
	#pragma omp parallel 
#endif
	{
		std::vector<arrayType<maxAlleles> > markerProbabilities(nFinerPoints);
		//These are inputs, used when we want to convert haplotype probabilities (n * n tables) into marker probabilities (these are much smaller tables typically 2*2)
		funnelHaplotypeToMarker<nFounders, maxAlleles> finerFunnelHaplotypeToMarker(finerFunnelHaplotypeProbabilities);
		funnelHaplotypeToMarker<nFounders, maxAlleles> normalFunnelHaplotypeToMarker(funnelHaplotypeProbabilities);

		intercrossingHaplotypeToMarker<nFounders, maxAlleles> finerIntercrossingHaplotypeToMarker(finerIntercrossingHaplotypeProbabilities);
		intercrossingHaplotypeToMarker<nFounders, maxAlleles> normalIntercrossingHaplotypeToMarker(intercrossingHaplotypeProbabilities);
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
						arrayType<maxAlleles>* markerProbabilitiesThisIntercrossing = &(thisMarkerPairData.perAIGenerationData(0, intercrossingGeneration-1, selfingCounter - minSelfing));
						if(thisMarkerPairData.allowableAI(intercrossingGeneration-1, selfingCounter - minSelfing))
						{
							normalIntercrossingHaplotypeToMarker.template convert<true>(markerProbabilitiesThisIntercrossing, intercrossingGeneration, firstMarkerPatternData, secondMarkerPatternData, selfingCounter);
						}
					}
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						arrayType<maxAlleles>* markerProbabilitiesThisFunnel = &(thisMarkerPairData.perFunnelData(0, funnelCounter, selfingCounter - minSelfing));
						memset(&markerProbabilitiesThisFunnel, 0, sizeof(arrayType<maxAlleles>));
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