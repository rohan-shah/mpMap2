#ifndef PROBABILITIES_8_HEADER_GUARD
#define PROBABILITIES_8_HEADER_GUARD
#include "probabilities.hpp"
template<> struct probabilityData<8>
{
public:
	/*See Karl Bromans paper on intermediate generations. This mask converts allele encodings (0 - 2) into indices into 
	the array of 4 different probabilities. In terms of indices, 
	Zero = homozygote, One = other homozygote, Two = hetrozygote
	In terms of values, see Table one of the paper. Zero = first equation of table, one = second equation, etc. Note that
	we combine equations 4 and 5 into a single state. 
	*/
	static const int intermediateProbabilitiesMask[64][64];
	/*This mask takes in the two alleles at a *single* location and returns a value encoding that genotype. */
	static const int intermediateAllelesMask[8][8];
	static const int infiniteMask[8][8];
private:
	probabilityData(){}
};
#endif
