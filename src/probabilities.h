#ifndef PROBABILITIES_HEADER_GUARD
#define PROBABILITIES_HEADER_GUARD
#include "matrices.h"
#include <cmath>
#include <stdexcept>
#include <limits>
#include <array>
#include "array2.h"
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
#include "compressedProbabilities.h"
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
template<int nFounders, bool infiniteSelfing, bool takeLogs> struct expandedGenotypeProbabilitiesPhased;
/*
 * First the infinite generations of selfing case
 */
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilitiesPhased<nFounders, true, takeLogs>
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
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilitiesPhased<nFounders, false, takeLogs>
{
public:
	static void noIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesNoIntercross<nFounders, false>(probabilities, r, selfingGenerations, nFunnels);
		memset(expandedProbabilities.values, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>::values));
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
		memset(expandedProbabilities.values, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>::values));
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
template<int nFounders, bool infiniteSelfing, bool takeLogs> struct expandedGenotypeProbabilitiesUnphased;
/*
 * First the infinite generations of selfing case. Phasing isn't an issue here, so this is the same as the unphased case. 
 */
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilitiesUnphased<nFounders, true, takeLogs>
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
 * Now the finite generations of selfing case. This case is *different* from the phased case above. 
 */
template<int nFounders, bool takeLogs> struct expandedGenotypeProbabilitiesUnphased<nFounders, false, takeLogs>
{
public:
	static void noIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesNoIntercross<nFounders, false>(probabilities, r, selfingGenerations, nFunnels);
		memset(expandedProbabilities.values, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>::values));
#ifdef INTERNAL_CHECKS
		double sum = 0;
#endif
		for(int marker1Allele1 = 0; marker1Allele1 < nFounders; marker1Allele1++)
		{
			for(int marker1Allele2 = 0; marker1Allele2 <= marker1Allele1; marker1Allele2++)
			{
				const int marker1Index1 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele1][marker1Allele2];
				const int marker1Index2 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele2][marker1Allele1];
				int multiple1 = 1;
				if(marker1Allele1 != marker1Allele2) multiple1 *= 2;
				for(int marker2Allele1 = 0; marker2Allele1 < nFounders; marker2Allele1++)
				{
					for(int marker2Allele2 = 0; marker2Allele2 <= marker2Allele1; marker2Allele2++)
					{
						int multiple2 = 1;
						if(marker2Allele1 != marker2Allele2) multiple2 *= 2;
						const int marker2Index1 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele1][marker2Allele2];
						const int marker2Index2 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele2][marker2Allele1];
						double value = 0.25 * (probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index1][marker2Index1]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index1][marker2Index2]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index2][marker2Index1]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index2][marker2Index2]]);
#ifdef INTERNAL_CHECKS
						sum += value * multiple1 * multiple2;
#endif
						if(takeLogs) value = log(value);
						expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2] = expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele2][marker2Allele1] = expandedProbabilities.values[marker1Allele2][marker1Allele1][marker2Allele1][marker2Allele2] = expandedProbabilities.values[marker1Allele2][marker1Allele1][marker2Allele2][marker2Allele1] = value;
					}
				}
			}
		}
#ifdef INTERNAL_CHECKS
		if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
	static void withIntercross(expandedProbabilitiesFiniteSelfing<nFounders>& expandedProbabilities, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
	{
		const int nDifferentProbs = compressedProbabilities<nFounders, false>::nDifferentProbs;
		std::array<double, nDifferentProbs> probabilities;
		genotypeProbabilitiesWithIntercross<nFounders, false>(probabilities, nAIGenerations, r, selfingGenerations, nFunnels);
		memset(expandedProbabilities.values, 0, sizeof(expandedProbabilitiesFiniteSelfing<nFounders>::values));
#ifdef INTERNAL_CHECKS
		double sum = 0;
#endif
		for(int marker1Allele1 = 0; marker1Allele1 < nFounders; marker1Allele1++)
		{
			for(int marker1Allele2 = 0; marker1Allele2 <= marker1Allele1; marker1Allele2++)
			{
				const int marker1Index1 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele1][marker1Allele2];
				const int marker1Index2 = probabilityData<nFounders>::intermediateAllelesMask[marker1Allele2][marker1Allele1];
				int multiple1 = 1;
				if(marker1Allele1 != marker1Allele2) multiple1 *= 2;
				for(int marker2Allele1 = 0; marker2Allele1 < nFounders; marker2Allele1++)
				{
					for(int marker2Allele2 = 0; marker2Allele2 <= marker2Allele1; marker2Allele2++)
					{
						int multiple2 = 1;
						if(marker2Allele1 != marker2Allele2) multiple2 *= 2;
						const int marker2Index1 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele1][marker2Allele2];
						const int marker2Index2 = probabilityData<nFounders>::intermediateAllelesMask[marker2Allele2][marker2Allele1];
						double value = 0.25 * (probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index1][marker2Index1]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index1][marker2Index2]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index2][marker2Index1]] + probabilities[probabilityData<nFounders>::intermediateProbabilitiesMask[marker1Index2][marker2Index2]]);
#ifdef INTERNAL_CHECKS
						sum += value * multiple1 * multiple2;
#endif
						if(takeLogs) value = log(value);
						expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele1][marker2Allele2] = expandedProbabilities.values[marker1Allele1][marker1Allele2][marker2Allele2][marker2Allele1] = expandedProbabilities.values[marker1Allele2][marker1Allele1][marker2Allele1][marker2Allele2] = expandedProbabilities.values[marker1Allele2][marker1Allele1][marker2Allele2][marker2Allele1] = value;
					}
				}
			}
		}
#ifdef INTERNAL_CHECKS
		if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Haplotype probabilities did not sum to 1");
#endif
	}
};
#endif
