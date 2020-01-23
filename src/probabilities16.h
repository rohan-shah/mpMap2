#ifndef PROBABILITIES_16_HEADER_GUARD
#define PROBABILITIES_16_HEADER_GUARD
#include "probabilities.h"
template<> struct probabilityData<16>
{
public:
	static const int intermediateProbabilitiesMask[256][256];
	static const int intermediateAllelesMask[16][16];
	static const int infiniteMask[16][16];
private:
	probabilityData(){}
};
#endif
