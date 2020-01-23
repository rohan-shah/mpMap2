#ifndef PROBABILITIES_2_HEADER_GUARD
#define PROBABILITIES_2_HEADER_GUARD
#include "probabilities.h"
template<> struct probabilityData<2>
{
public:
	static const int intermediateProbabilitiesMask[4][4];
	static const int intermediateAllelesMask[2][2];
	static const int infiniteMask[2][2];
private:
	probabilityData(){}
};
#endif
