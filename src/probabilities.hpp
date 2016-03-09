#ifndef PROBABILITIES_HEADER_GUARD
#define PROBABILITIES_HEADER_GUARD
#include <matrices.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>
const int nDifferentProbs = 100;
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
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesNoIntercross(double(&prob)[nDifferentProbs], double recombinationFraction, int selfingGenerations, std::size_t nFunnels, int& nValues);
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesWithIntercross(double(&prob)[nDifferentProbs], int nAIGenarations, double recombinationFraction, int selfingGenerations, std::size_t nFunnels, int& nValues);
//Templated function to work out the single locus probabilities. This is needed for imputation / the Viterbi algorithm
template<int nFounders, bool infiniteSelfing> void singleLocusGenotypeProbabilitiesNoIntercross(array2<nFounders>& data, int selfingGenerations, std::size_t nFunnels);
template<int nFounders, bool infiniteSelfing> void singleLocusGenotypeProbabilitiesWithIntercross(array2<nFounders>& data, int nAIGenarations, int selfingGenerations, std::size_t nFunnels);
template<int nFounders, bool infiniteSelfing, bool takeLogs> struct expandedGenotypeProbabilities;
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilities<nFounders, true, takeLogs>
{
public:
	static void noIntercross(array2<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		double probabilities[nDifferentProbs];
		int nValues;
		genotypeProbabilitiesNoIntercross<nFounders, true>(probabilities, r, selfingGenerations, nFunnels, nValues);
		if(takeLogs)
		{
			for(int i = 0; i < nValues; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				expandedProbabilities.values[i][j] = probabilities[probabilityData<nFounders>::infiniteMask[i][j]];
			}
		}
	}
	static void withIntercross(array2<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
	{
		double probabilities[nDifferentProbs];
		int nValues;
		genotypeProbabilitiesWithIntercross<nFounders, true>(probabilities, nAIGenerations, r, selfingGenerations, nFunnels, nValues);
		if(takeLogs)
		{
			for(int i = 0; i < nValues; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				expandedProbabilities.values[i][j] = probabilities[probabilityData<nFounders>::infiniteMask[i][j]];
			}
		}
	}
};
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilities<nFounders, false, takeLogs>
{
public:
	static void noIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		double probabilities[nDifferentProbs];
		int nValues;
		genotypeProbabilitiesNoIntercross<nFounders, false>(probabilities, r, selfingGenerations, nFunnels, nValues);
		memset(&expandedProbabilities, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>));
		if(takeLogs)
		{
			for(int i = 0; i < nValues; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
#ifndef NDEBUG
		double sum = 0;
#endif
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
#ifndef NDEBUG
						sum += expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2];
#endif
					}
				}
			}
		}
#ifndef NDEBUG
		if(!takeLogs && fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
	static void withIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
	{	
		double probabilities[nDifferentProbs];
		int nValues;
		genotypeProbabilitiesWithIntercross<nFounders, false>(probabilities, nAIGenerations, r, selfingGenerations, nFunnels, nValues);
		memset(&expandedProbabilities, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>));
		if(takeLogs)
		{
			for(int i = 0; i < nValues; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
#ifndef NDEBUG
		double sum = 0;
#endif
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
#ifndef NDEBUG
						sum += expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2];
#endif
					}
				}
			}
		}
#ifndef NDEBUG
		if(!takeLogs && fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
};
#endif
