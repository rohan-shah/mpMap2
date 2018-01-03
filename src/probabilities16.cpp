#include "probabilities16.h"
#include <cmath>
#include <stdexcept>
const int probabilityData<16>::intermediateProbabilitiesMask[][256] = 
#include "probabilities16IntermediateProbabilitiesMask.h"
;
//Matrix used to encode a pair of founders into a genotype.
const int probabilityData<16>::intermediateAllelesMask[][16] = 
  {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {16, 17, 18, 
  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}, {32, 33, 34, 
  35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47}, {48, 49, 50, 
  51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}, {64, 65, 66, 
  67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79}, {80, 81, 82, 
  83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95}, {96, 97, 98, 
  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 
  111}, {112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 
  124, 125, 126, 127}, {128, 129, 130, 131, 132, 133, 134, 135, 136, 
  137, 138, 139, 140, 141, 142, 143}, {144, 145, 146, 147, 148, 149, 
  150, 151, 152, 153, 154, 155, 156, 157, 158, 159}, {160, 161, 162, 
  163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 
  175}, {176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 
  188, 189, 190, 191}, {192, 193, 194, 195, 196, 197, 198, 199, 200, 
  201, 202, 203, 204, 205, 206, 207}, {208, 209, 210, 211, 212, 213, 
  214, 215, 216, 217, 218, 219, 220, 221, 222, 223}, {224, 225, 226, 
  227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 
  239}, {240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 
  252, 253, 254, 255}};
//In the case of infinite generations of selfing, this table is used to turn the computed probabilities (there are four unique values) into a 16 x 16 probability matrix (containing only those three unique values)
const int probabilityData<16>::infiniteMask[][16] = 
		{{0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		{1, 0, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		{2, 2, 0, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		{2, 2, 1, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		{3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3},
		{3, 3, 3, 3, 1, 0, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3},
		{3, 3, 3, 3, 2, 2, 0, 1, 3, 3, 3, 3, 3, 3, 3, 3},
		{3, 3, 3, 3, 2, 2, 1, 0, 3, 3, 3, 3, 3, 3, 3, 3},
		{3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3},
		{3, 3, 3, 3, 3, 3, 3, 3, 1, 0, 2, 2, 3, 3, 3, 3},
		{3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 1, 3, 3, 3, 3},
		{3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 0, 3, 3, 3, 3},
		{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2},
		{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 0, 2, 2},
		{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 1},
		{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 0}};
template<> void genotypeProbabilitiesNoIntercross<16, true>(std::array<double, 4>& prob, double r, int, std::size_t)
{
	prob[0] = (1-r)*(1-r)*(1-r)/(16*(1 + 2*r));
	prob[1] = (1-r)*(1-r)*r/(16+32*r);
	prob[2] = (1 - r)*r/(32+64*r);
	prob[3] = r/(64*(1 + 2 * r));
}
template<> void genotypeProbabilitiesNoIntercross<16, false>(std::array<double, 95>& prob, double r, int selfingGenerations, std::size_t nFunnels)
{
	double onePlus2R = 1 + 2 *r;
	double rSquared = r*r;
	double oneMinusRSquared = (1-r)*(1-r);
	double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
	double powD1 = std::pow(1 + 2*(-1 + r)*r, selfingGenerations);
	double oneMinusRCubed = (1-r)*oneMinusRSquared;
	double oneMinusRPow4 = oneMinusRSquared*oneMinusRSquared;
	double oneMinusRPow5 = oneMinusRCubed*oneMinusRSquared;
	double oneMinusRPow6 = oneMinusRCubed*oneMinusRCubed;
	double oneMinusR = 1 - r;
	double pow2 = std::pow(2, selfingGenerations);
    
	prob[0] = -(oneMinusRCubed*(2 - 2*pow2 - onePlus2R*powD1 + powOneMinus2R + 4*r - 2*powOneMinus2R*r))/(32*onePlus2R*pow2);
	prob[1] = 0;
	prob[2] = 0;
	prob[3] = 0;
	prob[4] = -(oneMinusRCubed*(-1 + powD1))/(256*pow2);
	prob[5] = (oneMinusRSquared*r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(32*onePlus2R*pow2);
	prob[6] = 0;
	prob[7] = 0;
	prob[8] = -(oneMinusRSquared*(-1 + powD1)*r)/(256*pow2);
	prob[9] = (oneMinusR*r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(64*onePlus2R*pow2);
	prob[10] = 0;
	prob[11] = 0;
	prob[12] = -(oneMinusR*(-1 + powD1)*r)/(512*pow2);
	prob[13] = (r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(128*onePlus2R*pow2);
	prob[14] = 0;
	prob[15] = 0;
	prob[16] = -((-1 + powD1)*r)/(1024*pow2);
	prob[17] = (-2 + onePlus2R*powD1 + powOneMinus2R - 4*r + 4*pow2*r - 2*powOneMinus2R*r)/(256*onePlus2R*pow2);
	prob[18] = 0;
	prob[19] = 0;
	prob[20] = 0;
	prob[21] = 0;
	prob[22] = 0;
	prob[23] = 0;
	prob[24] = 0;
	prob[25] = 0;
	prob[26] = 0;
	prob[27] = 0;
	prob[28] = 0;
	prob[29] = 0;
	prob[30] = 0;
	prob[31] = 0;
	prob[32] = 0;
	prob[33] = 0;
	prob[34] = 0;
	prob[35] = 0;
	prob[36] = 0;
	prob[37] = 0;
	prob[38] = 0;
	prob[39] = 0;
	prob[40] = 0;
	prob[41] = 0;
	prob[42] = 0;
	prob[43] = 0;
	prob[44] = 0;
	prob[45] = 0;
	prob[46] = 0;
	prob[47] = 0;
	prob[48] = 0;
	prob[49] = 0;
	prob[50] = 0;
	prob[51] = 0;
	prob[52] = 0;
	prob[53] = 0;
	prob[54] = 0;
	prob[55] = 0;
	prob[56] = 0;
	prob[57] = 0;
	prob[58] = 0;
	prob[59] = 0;
	prob[60] = 0;
	prob[61] = 0;
	prob[62] = 0;
	prob[63] = 0;
	prob[64] = 0;
	prob[65] = 0;
	prob[66] = 0;
	prob[67] = 0;
	prob[68] = 0;
	prob[69] = 0;
	prob[70] = 0;
	prob[71] = 0;
	prob[72] = 0;
	prob[73] = 0;
	prob[74] = 0;
	prob[75] = (oneMinusRPow6*(powD1 + powOneMinus2R))/(256*pow2);
	prob[76] = (oneMinusRPow5*(powD1 + powOneMinus2R)*r)/(256*pow2);
	prob[77] = (oneMinusRPow4*(powD1 + powOneMinus2R)*r)/(512*pow2);
	prob[78] = (oneMinusRCubed*(powD1 + powOneMinus2R)*r)/(1024*pow2);
	prob[79] = (oneMinusRPow4*(powD1 + powOneMinus2R)*rSquared)/(256*pow2);
	prob[80] = (oneMinusRCubed*(powD1 + powOneMinus2R)*rSquared)/(512*pow2);
	prob[81] = (oneMinusRSquared*(powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
	prob[82] = (oneMinusRSquared*(powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
	prob[83] = (oneMinusR*(powD1 + powOneMinus2R)*rSquared)/(2048*pow2);
	prob[84] = ((powD1 + powOneMinus2R)*rSquared)/(4096*pow2);
	prob[85] = -(oneMinusRPow6*(-powD1 + powOneMinus2R))/(256*pow2);
	prob[86] = -(oneMinusRPow5*(-powD1 + powOneMinus2R)*r)/(256*pow2);
	prob[87] = -(oneMinusRPow4*(-powD1 + powOneMinus2R)*r)/(512*pow2);
	prob[88] = -(oneMinusRCubed*(-powD1 + powOneMinus2R)*r)/(1024*pow2);
	prob[89] = -(oneMinusRPow4*(-powD1 + powOneMinus2R)*rSquared)/(256*pow2);
	prob[90] = -(oneMinusRCubed*(-powD1 + powOneMinus2R)*rSquared)/(512*pow2);
	prob[91] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
	prob[92] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
	prob[93] = -(oneMinusR*(-powD1 + powOneMinus2R)*rSquared)/(2048*pow2);
	prob[94] = ((powD1 - powOneMinus2R)*rSquared)/(4096*pow2);
}
template<> void genotypeProbabilitiesWithIntercross<16, true>(std::array<double, 4>& prob, int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)*(1-r)*(1-r)/16 + (2 * r + 1 - powOneMinusR)/(16*16))/(1 + 2 * r);
	prob[1] = prob[2] = prob[3] = (1 - 16*prob[0])/(16*15);
}
template<> void genotypeProbabilitiesWithIntercross<16, false>(std::array<double, 95>& prob, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		throw std::runtime_error("The case of a single funnel with intercrossing is not yet handled");
	}
  //See mathematica code for rewritting rules
	double oneMinusR = (1 - r);
	double onePlus2R = 1 + 2*r;
	double oneMinusRSquared = oneMinusR*oneMinusR;
	double twoRMinus3 = 2 * r - 3;
	double twoRMinus3Squared = twoRMinus3*twoRMinus3;
	double powOneMinusR1 = std::pow(oneMinusR, 1+nAIGenerations);
	double powOneMinusR2 = std::pow(oneMinusR, nAIGenerations - 1);
	double powOneMinusR3 = std::pow(oneMinusR, nAIGenerations - 2);
	double oneMinusRPow1 = std::pow(oneMinusR, 2*nAIGenerations);
	double pow2 = std::pow(2, selfingGenerations);
	double twoRMinus1 = 2*r - 1;
	double complexPower1 = std::pow(1 + 2 *(-1 + r)* r, selfingGenerations);
	double quadratic1 = 5 + 4 * (-2 + r) * r;
	double quadratic1Squared = quadratic1*quadratic1;
	double twoRMinus1Squared = twoRMinus1*twoRMinus1;
	double twoRMinus1PowD = std::pow(-twoRMinus1, selfingGenerations);
	double rMinus1Pow4 = std::pow(r - 1, 4);
	double complexPart1 = 1.0/256.0 + powOneMinusR2 * (-1.0/256.0 + (1.0/16.0) *rMinus1Pow4);
	double complexPart2 = -7 + std::pow(2.0, 3.0 + selfingGenerations) + 4*twoRMinus1*twoRMinus1PowD - 14*r;
	double complexPart3 = 1.0 - 16.0*complexPart1;
	double complexPart3Squared = complexPart3*complexPart3;
	double complexPart1Squared = complexPart1*complexPart1;
	double complexPart4 = complexPower1 + 8*powOneMinusR2*twoRMinus1PowD;
	double complexPart5 = 105*complexPower1 - 16*(-7*complexPower1 + 8*twoRMinus1PowD)*(-2 + r)*r*(2 + (-2 + r)*r);
	double complexPart6 = complexPart3*(15*complexPart4 + powOneMinusR2*(105*complexPower1 + 16*(-2 + r)*r*(2 + (-2 + r)*r)*(7*complexPower1 + 8*twoRMinus1PowD)));
	double complexPart7 = complexPart3*(complexPart5*powOneMinusR2 + 15*(complexPower1 - 8*powOneMinusR2*twoRMinus1PowD));
	double complexPart8 = 15*(-1 + r) + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3;
	double complexPart8Squared = complexPart8*complexPart8;
	double complexPart9 = complexPart8*(-240*oneMinusR + complexPower1*(225*oneMinusR + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3));
	double complexPart10 = 240*oneMinusRSquared + 112*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 - complexPower1*(225*oneMinusRSquared + 98*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared * twoRMinus3Squared);
	double complexPart11 = 57600*complexPart1Squared*(complexPower1 - twoRMinus1PowD) + complexPart3Squared*(complexPower1 + twoRMinus1PowD);
	double complexPart12 = complexPart3Squared*(complexPower1 - twoRMinus1PowD) + 57600*complexPart1Squared*(complexPower1 + twoRMinus1PowD);
	double complexPart13 = std::pow(0.5 - oneMinusR*r, selfingGenerations);
	double complexPart14 = -15*(-15 + std::pow(2.0, 3 + selfingGenerations))*oneMinusR *onePlus2R + complexPart2*oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3;
	double complexPart15 = 32*((-15 + std::pow(2.0, 3 + selfingGenerations))*oneMinusRSquared*onePlus2R + complexPart2*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3) + complexPower1*onePlus2R*(225*oneMinusRSquared + 98*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared * twoRMinus3Squared);
	double complexPart16 = -480*complexPart14*oneMinusR + complexPower1*onePlus2R*(50625*oneMinusRSquared - 1470*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared * twoRMinus3Squared);
	double complexPart17 = 240 - 225*complexPower1 + 1680*powOneMinusR2 + powOneMinusR3*(1792*oneMinusR*(-2 + r)*r*(2 + (-2 + r)*r) - complexPower1*quadratic1*twoRMinus1*twoRMinus3*(98*oneMinusR + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3));
	double complexPart18 = std::pow(4 - 8*oneMinusR*r, selfingGenerations);
	double complexPart19 = std::pow(2.0, -3*(1 + selfingGenerations));

	prob[0] = complexPart15/(4096*oneMinusRSquared*onePlus2R*pow2);
	prob[1] = complexPart17/(15360*pow2);
	prob[2] = complexPart10/(7680*oneMinusRSquared*pow2);
	prob[3] = complexPart10/(3840*oneMinusRSquared*pow2);
	prob[4] = complexPart10/(1920*oneMinusRSquared*pow2);
	prob[5] = complexPart16/(921600*oneMinusRSquared*onePlus2R*pow2);
	prob[6] = complexPart9/(115200*oneMinusRSquared*pow2);
	prob[7] = complexPart9/(57600*oneMinusRSquared*pow2);
	prob[8] = complexPart9/(28800*oneMinusRSquared*pow2);
	prob[9] = complexPart16/(460800*oneMinusRSquared*onePlus2R*pow2);
	prob[10] = complexPart9/(230400*oneMinusRSquared*pow2);
	prob[11] = complexPart9/(28800*oneMinusRSquared*pow2);
	prob[12] = complexPart9/(14400*oneMinusRSquared*pow2);
	prob[13] = complexPart16/(230400*oneMinusRSquared*onePlus2R*pow2);
	prob[14] = complexPart9/(115200*oneMinusRSquared*pow2);
	prob[15] = complexPart9/(57600*oneMinusRSquared*pow2);
	prob[16] = complexPart9/(7200*oneMinusRSquared*pow2);
	prob[17] = complexPart16/(115200*oneMinusRSquared*onePlus2R*pow2);
	prob[18] = complexPart9/(57600*oneMinusRSquared*pow2);
	prob[19] = complexPart9/(28800*oneMinusRSquared*pow2);
	prob[20] = complexPart9/(14400*oneMinusRSquared*pow2);
	prob[21] = complexPart12/(7200*pow2);
	prob[22] = complexPart6/(7200*pow2);
	prob[23] = complexPart6/(3600*pow2);
	prob[24] = complexPart6/(1800*pow2);
	prob[25] = complexPart11/(7200*pow2);
	prob[26] = complexPart7/(7200*pow2);
	prob[27] = complexPart7/(3600*pow2);
	prob[28] = complexPart7/(1800*pow2);
	prob[29] = (complexPart18*complexPart19*complexPart3Squared)/225;
	prob[30] = (2*complexPart3Squared*complexPower1)/(225*pow2);
	prob[31] = (4*complexPart3Squared*complexPower1)/(225*pow2);
	prob[32] = (complexPart3Squared*complexPower1)/(900*pow2);
	prob[33] = (complexPart13*complexPart3Squared)/225;
	prob[34] = (8*complexPart3Squared*complexPower1)/(225*pow2);
	prob[35] = (complexPart3Squared*complexPower1)/(450*pow2);
	prob[36] = (2*complexPart3Squared*complexPower1)/(225*pow2);
	prob[37] = (4*complexPart3Squared*complexPower1)/(225*pow2);
	prob[38] = complexPart12/(3600*pow2);
	prob[39] = complexPart6/(14400*pow2);
	prob[40] = complexPart6/(1800*pow2);
	prob[41] = complexPart6/(900*pow2);
	prob[42] = (complexPart8Squared*complexPower1)/(460800*oneMinusRSquared*pow2);
	prob[43] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[44] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
	prob[45] = complexPart11/(3600*pow2);
	prob[46] = complexPart7/(14400*pow2);
	prob[47] = complexPart7/(1800*pow2);
	prob[48] = complexPart7/(900*pow2);
	prob[49] = (complexPart8Squared*complexPower1)/(460800*oneMinusRSquared*pow2);
	prob[50] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[51] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
	prob[52] = (complexPart13*complexPart3Squared)/225;
	prob[53] = (16*complexPart3Squared*complexPower1)/(225*pow2);
	prob[54] = (2*complexPart3Squared*complexPower1)/(225*pow2);
	prob[55] = (8*complexPart3Squared*complexPower1)/(225*pow2);
	prob[56] = complexPart12/(1800*pow2);
	prob[57] = complexPart6/(7200*pow2);
	prob[58] = complexPart6/(3600*pow2);
	prob[59] = complexPart6/(450*pow2);
	prob[60] = (complexPart8Squared*complexPower1)/(230400*oneMinusRSquared*pow2);
	prob[61] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
	prob[62] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
	prob[63] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
	prob[64] = (complexPart8Squared*complexPower1)/(3600*oneMinusRSquared*pow2);
	prob[65] = complexPart11/(1800*pow2);
	prob[66] = complexPart7/(7200*pow2);
	prob[67] = complexPart7/(3600*pow2);
	prob[68] = complexPart7/(450*pow2);
	prob[69] = (complexPart8Squared*complexPower1)/(230400*oneMinusRSquared*pow2);
	prob[70] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
	prob[71] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
	prob[72] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
	prob[73] = (complexPart8Squared*complexPower1)/(3600*oneMinusRSquared*pow2);
	prob[74] = (8*complexPart3Squared*complexPower1)/(225*pow2);
	prob[75] = complexPart12/(900*pow2);
	prob[76] = complexPart6/(3600*pow2);
	prob[77] = complexPart6/(1800*pow2);
	prob[78] = complexPart6/(900*pow2);
	prob[79] = (complexPart8Squared*complexPower1)/(115200*oneMinusRSquared*pow2);
	prob[80] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[81] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
	prob[82] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[83] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
	prob[84] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
	prob[85] = complexPart11/(900*pow2);
	prob[86] = complexPart7/(3600*pow2);
	prob[87] = complexPart7/(1800*pow2);
	prob[88] = complexPart7/(900*pow2);
	prob[89] = (complexPart8Squared*complexPower1)/(115200*oneMinusRSquared*pow2);
	prob[90] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[91] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
	prob[92] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
	prob[93] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
	prob[94] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
#ifdef INTERNAL_CHECKS
	double sum = 0;
	for(int i = 0; i < 95; i++) sum += prob[i];
	if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Internal error");
#endif
	//This is because we combined some states (see mathematica code)
	prob[0] /= 16;
	prob[1] /= 64;
	prob[2] /= 128;
	prob[3] /= 256;
	prob[4] /= 512;
	prob[5] /= 16;
	prob[6] /= 128;
	prob[7] /= 256;
	prob[8] /= 512;
	prob[9] /= 32;
	prob[10] /= 64;
	prob[11] /= 512;
	prob[12] /= 1024;
	prob[13] /= 64;
	prob[14] /= 128;
	prob[15] /= 256;
	prob[16] /= 2048;
	prob[17] /= 128;
	prob[18] /= 256;
	prob[19] /= 512;
	prob[20] /= 1024;
	prob[21] /= 16;
	prob[22] /= 128;
	prob[23] /= 256;
	prob[24] /= 512;
	prob[25] /= 16;
	prob[26] /= 128;
	prob[27] /= 256;
	prob[28] /= 512;
	prob[29] /= 32;
	prob[30] /= 512;
	prob[31] /= 1024;
	prob[32] /= 64;
	prob[33] /= 256;
	prob[34] /= 2048;
	prob[35] /= 128;
	prob[36] /= 512;
	prob[37] /= 1024;
	prob[38] /= 32;
	prob[39] /= 64;
	prob[40] /= 512;
	prob[41] /= 1024;
	prob[42] /= 32;
	prob[43] /= 512;
	prob[44] /= 1024;
	prob[45] /= 32;
	prob[46] /= 64;
	prob[47] /= 512;
	prob[48] /= 1024;
	prob[49] /= 32;
	prob[50] /= 512;
	prob[51] /= 1024;
	prob[52] /= 256;
	prob[53] /= 4096;
	prob[54] /= 512;
	prob[55] /= 2048;
	prob[56] /= 64;
	prob[57] /= 128;
	prob[58] /= 256;
	prob[59] /= 2048;
	prob[60] /= 64;
	prob[61] /= 256;
	prob[62] /= 2048;
	prob[63] /= 256;
	prob[64] /= 4096;
	prob[65] /= 64;
	prob[66] /= 128;
	prob[67] /= 256;
	prob[68] /= 2048;
	prob[69] /= 64;
	prob[70] /= 256;
	prob[71] /= 2048;
	prob[72] /= 256;
	prob[73] /= 4096;
	prob[74] /= 2048;
	prob[75] /= 128;
	prob[76] /= 256;
	prob[77] /= 512;
	prob[78] /= 1024;
	prob[79] /= 128;
	prob[80] /= 512;
	prob[81] /= 1024;
	prob[82] /= 512;
	prob[83] /= 2048;
	prob[84] /= 2048;
	prob[85] /= 128;
	prob[86] /= 256;
	prob[87] /= 512;
	prob[88] /= 1024;
	prob[89] /= 128;
	prob[90] /= 512;
	prob[91] /= 1024;
	prob[92] /= 512;
	prob[93] /= 2048;
	prob[94] /= 2048;
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<16, false>(array2<16>&data, int selfingGenerations, std::size_t nFunnels)
{
	throw std::runtime_error("Single locus probabilities not implemented");
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<16, true>(array2<16>&data, int selfingGenerations, std::size_t nFunnels)
{
	throw std::runtime_error("Single locus probabilities not implemented");
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<16, false>(array2<16>& data, int selfingGenerations, std::size_t nFunnels)
{
	throw std::runtime_error("Single locus probabilities not implemented");
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<16, true>(array2<16>& data, int selfingGenerations, std::size_t nFunnels)
{
	throw std::runtime_error("Single locus probabilities not implemented");
}
