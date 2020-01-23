#ifndef PROBABILITIES_8_HEADER_GUARD
#define PROBABILITIES_8_HEADER_GUARD
#include "probabilities.h"
template<> struct probabilityData<8>
{
public:
	static const int intermediateProbabilitiesMask[64][64];
	static const int intermediateAllelesMask[8][8];
	static const int infiniteMask[8][8];
private:
	probabilityData(){}
};
#endif
