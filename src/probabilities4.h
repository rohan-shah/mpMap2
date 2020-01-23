#ifndef PROBABILITIES_4_HEADER_GUARD
#define PROBABILITIES_4_HEADER_GUARD
#include "probabilities.h"
template<> struct probabilityData<4>
{
public:
	static const int intermediateProbabilitiesMask[16][16];
	static const int intermediateAllelesMask[4][4];
	static const int infiniteMask[4][4];
private:
	probabilityData(){}
};
#endif
