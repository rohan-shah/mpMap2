#ifndef PROBABILITIES_HEADER_GUARD
#define PROBABILITIES_HEADER_GUARD
#include <matrices.hpp>

const int nDifferentProbs = 10;
template<int nFounders> struct probabilityData;
template<int nFounders> struct expandedProbabilitiesFiniteSelfing
{
public:
	expandedProbabilitiesFiniteSelfing()
	{}
	double values[nFounders][nFounders][nFounders][nFounders];
};
template<int nFounders, bool infiniteSelfing> struct expandedProbabilities;
template<int nFounders> struct expandedProbabilities<nFounders, true>
{
	typedef array2<nFounders> type;
};
template<int nFounders> struct expandedProbabilities<nFounders, false>
{
	typedef expandedProbabilitiesFiniteSelfing<nFounders> type;
};

//Templated function to work out the two-point probabilities with the given recombination fraction (and number of AI generations). Templating allows the number of founders to be a compile-time constant
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesNoIntercross(double (&prob)[nDifferentProbs], double recombinationFraction, int selfingGenerations);
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesWithIntercross(double (&prob)[nDifferentProbs], int nAIGenarations, double recombinationFraction, int selfingGenerations);
template<int nFounders, bool infiniteSelfing> struct expandedGenotypeProbabilities;
template<int nFounders> struct expandedGenotypeProbabilities<nFounders, true>
{
public:
	static void noIntercross(array2<nFounders>& expandedProbabilities, double r, int selfingGenerations)
	{
		double probabilities[nDifferentProbs];
		genotypeProbabilitiesNoIntercross<nFounders, true>(probabilities, r, selfingGenerations);
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				expandedProbabilities.values[i][j] = probabilities[probabilityData<nFounders>::infiniteMask[i][j]];
			}
		}
	}
	static void withIntercross(array2<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations)
	{
		double probabilities[nDifferentProbs];
		genotypeProbabilitiesWithIntercross<nFounders, true>(probabilities, nAIGenerations, r, selfingGenerations);
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				expandedProbabilities.values[i][j] = probabilities[probabilityData<nFounders>::infiniteMask[i][j]];
			}
		}
	}
};
template<int nFounders> struct expandedGenotypeProbabilities<nFounders, false>
{
public:
	static void noIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, double r, int selfingGenerations)
	{
		double probabilities[nDifferentProbs];
		genotypeProbabilitiesNoIntercross<nFounders, false>(probabilities, r, selfingGenerations);
		for(int marker1Allele1 = 0; marker1Allele1 < nFounders; marker1Allele1++)
		{
			for(int marker1Allele2 = 0; marker1Allele2 < nFounders; marker1Allele2++)
			{
				for(int marker2Allele1 = 0; marker2Allele1 < nFounders; marker2Allele1++)
				{
					for(int marker2Allele2 = 0; marker2Allele2 < nFounders; marker2Allele2++)
					{
						const int index1 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele1][marker1Allele2];
						const int index2 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele1][marker2Allele2];
						expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2] = probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[index1][index2]];
					}
				}
			}
		}
	}
	static void withIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations)
	{	
		double probabilities[nDifferentProbs];
		genotypeProbabilitiesWithIntercross<nFounders, false>(probabilities, nAIGenerations, r, selfingGenerations);
		for(int marker1Allele1 = 0; marker1Allele1 < nFounders; marker1Allele1++)
		{
			for(int marker1Allele2 = 0; marker1Allele2 < nFounders; marker1Allele2++)
			{
				for(int marker2Allele1 = 0; marker2Allele1 < nFounders; marker2Allele1++)
				{
					for(int marker2Allele2 = 0; marker2Allele2 < nFounders; marker2Allele2++)
					{
						const int index1 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele1][marker1Allele2];
						const int index2 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele1][marker2Allele2];
						expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2] = probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[index1][index2]];
					}
				}
			}
		}
	}
};
//template<int nFounders, bool infiniteSelfing> void expandedGenotypeProbabilitiesNoIntercross(typename expandedProbabilities<nFounders, infiniteSelfing>::type& expandedProbabilities, double r, int selfingGenerations);
//template<int nFounders, bool infiniteSelfing> void expandedGenotypeProbabilitiesWithIntercross(typename expandedProbabilities<nFounders, infiniteSelfing>::type& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations);
#endif