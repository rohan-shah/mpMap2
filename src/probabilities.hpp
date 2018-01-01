#ifndef PROBABILITIES_HEADER_GUARD
#define PROBABILITIES_HEADER_GUARD
#include "matrices.hpp"
#include <cmath>
#include <stdexcept>
#include <limits>
#include <array>
#include "array2.hpp"
/*
 * Struct that will contain arrays relevant for probability calculations
 */
template<int nFounders> struct probabilityData;
/*
 * The type for the expanded probability data, in the case of infinite selfing
 */
template<int nFounders> struct expandedProbabilitiesFiniteSelfing
{
public:
	expandedProbabilitiesFiniteSelfing()
	{}
	double values[nFounders][nFounders][nFounders][nFounders];
};
/*
 * Type with a typedef, giving the type for the expanded probability data, both infinite generations of selfing and finite generations of selfing
 */
template<int nFounders, bool infiniteSelfing> struct expandedProbabilities;
/*
 * In the case of infinite selfing it's just a two-dimensional array
 */
template<int nFounders> struct expandedProbabilities<nFounders, true>
{
	typedef array2<nFounders> type;
};
/*
 * In the finite selfing case we use the struct specially define above
 */
template<int nFounders> struct expandedProbabilities<nFounders, false>
{
	typedef expandedProbabilitiesFiniteSelfing<nFounders> type;
};
#include "compressedProbabilities.hpp"
/*
 * Templated function to work out the two-point probabilities with the given recombination fraction (and number of AI generations). Templating allows the number of founders to be a compile-time constant
 */
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesNoIntercross(std::array<double, compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs>& prob, double recombinationFraction, int selfingGenerations, std::size_t nFunnels);
template<int nFounders, bool infiniteSelfing> void genotypeProbabilitiesWithIntercross(std::array<double, compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs>& prob, int nAIGenerations, double recombinationFraction, int selfingGenerations, std::size_t nFunnels);
/*
 * Templated function to work out the single locus probabilities. This is needed for imputation / the Viterbi algorithm
 */
template<int nFounders, bool infiniteSelfing> void singleLocusGenotypeProbabilitiesNoIntercross(array2<nFounders>& data, int selfingGenerations, std::size_t nFunnels);
template<int nFounders, bool infiniteSelfing> void singleLocusGenotypeProbabilitiesWithIntercross(array2<nFounders>& data, int selfingGenerations, std::size_t nFunnels);
/*
 * This struct has a static function that computes the compressed probabilities and immediately exands them to the right type
 */ 
template<int nFounders, bool infiniteSelfing, bool takeLogs> struct expandedGenotypeProbabilities;
/*
 * First the infinite generations of selfing case
 */
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilities<nFounders, true, takeLogs>
{
public:
	static void noIntercross(array2<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, true>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesNoIntercross<nFounders, true>(probabilities, r, selfingGenerations, nFunnels);
		if(takeLogs)
		{
			for(int i = 0; i < nDifferentProbs; i++)
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
		const int nDifferentProbs = compressedProbabilities<nFounders, true>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesWithIntercross<nFounders, true>(probabilities, nAIGenerations, r, selfingGenerations, nFunnels);
		if(takeLogs)
		{
			for(int i = 0; i < nDifferentProbs; i++)
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
/*
 * Now the finite generations of selfing case
 */
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilities<nFounders, false, takeLogs>
{
public:
	static void noIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesNoIntercross<nFounders, false>(probabilities, r, selfingGenerations, nFunnels);
		memset(&expandedProbabilities, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>));
		if(takeLogs)
		{
			for(int i = 0; i < nDifferentProbs; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
#ifdef INTERNAL_CHECKS
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
#ifdef INTERNAL_CHECKS
						sum += expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2];
#endif
					}
				}
			}
		}
#ifdef INTERNAL_CHECKS
		if(!takeLogs && fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
	static void withIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesWithIntercross<nFounders, false>(probabilities, nAIGenerations, r, selfingGenerations, nFunnels);
		memset(&expandedProbabilities, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>));
		if(takeLogs)
		{
			for(int i = 0; i < nDifferentProbs; i++)
			{
				if(probabilities[i] == 0) probabilities[i] = -std::numeric_limits<double>::infinity();
				else probabilities[i] = log(probabilities[i]);
			}
		}
#ifdef INTERNAL_CHECKS
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
#ifdef INTERNAL_CHECKS
						sum += expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2];
#endif
					}
				}
			}
		}
#ifdef INTERNAL_CHECKS
		if(!takeLogs && fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
};
#endif
