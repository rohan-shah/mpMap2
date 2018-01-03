#include "probabilities4.h"
#include <cmath>
#include <stdexcept>
//Matrix used to turn a pair of genotypes into a probability. 
const int probabilityData<4>::intermediateProbabilitiesMask[][16] =
{ { 0, 1, 2, 2, 1, 3, 4, 4, 2, 4, 5, 6, 2, 4, 6, 5 }, { 1, 7, 8, 8, 9, 1,
10, 10, 10, 8, 6, 11, 10, 8, 11, 6 }, { 2, 8, 12, 13, 10, 4, 13, 14,
15, 16, 2, 10, 16, 17, 8, 4 }, { 2, 8, 13, 12, 10, 4, 14, 13, 16, 17,
4, 8, 15, 16, 10, 2 }, { 1, 9, 10, 10, 7, 1, 8, 8, 8, 10, 6, 11, 8,
10, 11, 6 }, { 3, 1, 4, 4, 1, 0, 2, 2, 4, 2, 5, 6, 4, 2, 6, 5 }, { 4,
10, 13, 14, 8, 2, 12, 13, 16, 15, 2, 10, 17, 16, 8, 4 }, { 4, 10, 14,
13, 8, 2, 13, 12, 17, 16, 4, 8, 16, 15, 10, 2 }, { 2, 10, 15, 16, 8,
4, 16, 17, 12, 13, 2, 8, 13, 14, 10, 4 }, { 4, 8, 16, 17, 10, 2, 15,
16, 13, 12, 2, 8, 14, 13, 10, 4 }, { 5, 6, 2, 4, 6, 5, 2, 4, 2, 2, 0,
1, 4, 4, 1, 3 }, { 6, 11, 10, 8, 11, 6, 10, 8, 8, 8, 1, 7, 10, 10, 9,
1 }, { 2, 10, 16, 15, 8, 4, 17, 16, 13, 14, 4, 10, 12, 13, 8, 2 }, { 4,
8, 17, 16, 10, 2, 16, 15, 14, 13, 4, 10, 13, 12, 8, 2 }, { 6, 11, 8,
10, 11, 6, 8, 10, 10, 10, 1, 9, 8, 8, 7, 1 }, { 5, 6, 4, 2, 6, 5, 4,
2, 4, 4, 3, 1, 2, 2, 1, 0 } };
//Matrix used to encode a pair of founders into a genotype.
const int probabilityData<4>::intermediateAllelesMask[][4] = 
		{
			{0,  1,  2,  3},
			{4,  5,  6,  7},
			{8,  9,  10, 11},
			{12, 13, 14, 15}
		};
//In the case of infinite generations of selfing, this table is used to turn the computed probabilities (there are three computed values, two unique) into a 4 x 4 probability matrix (containing only those three values)
const int probabilityData<4>::infiniteMask[][4] = 
		{
			{0, 1, 2, 2},
			{1, 0, 2, 2},
			{2, 2, 0, 1},
			{2, 2, 1, 0}
		};
template<> void genotypeProbabilitiesNoIntercross<4, true>(std::array<double, 3>& prob, double r, int, std::size_t)
{
	prob[0] = (1-r)/(4*(1 + 2*r));
	prob[1] = prob[2] = r/(4*(1 + 2 * r));
}
template<> void genotypeProbabilitiesNoIntercross<4, false>(std::array<double, 18>& prob, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		double oneMinusR = 1 - r;
		double oneMinusRSquared = oneMinusR*oneMinusR;
		double onePlus2R = 1 + 2 * r;
		double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
		double powD1 = std::pow(1 + 2 * (-1 + r)*r, selfingGenerations);
		double pow2 = std::pow(2, selfingGenerations);
		double rSquared = r*r;
		prob[0] = (oneMinusR*(-2 + 2 * pow2 + onePlus2R*powD1 - powOneMinus2R - 4 * r + 2 * powOneMinus2R*r)) / (8 * onePlus2R*pow2);
		prob[1] = 0;
		prob[2] = -(oneMinusR*(-1 + powD1)) / (16 * pow2);
		prob[3] = (r*(-2 + 2 * pow2 + onePlus2R*powD1 - powOneMinus2R - 4 * r + 2 * powOneMinus2R*r)) / (8 * onePlus2R*pow2);
		prob[4] = -((-1 + powD1)*r) / (16 * pow2);
		prob[5] = (-2 + onePlus2R*powD1 + powOneMinus2R - 4 * r + 4 * pow2*r - 2 * powOneMinus2R*r) / (16 * onePlus2R*pow2);
		prob[6] = 0;
		prob[7] = 0;
		prob[8] = 0;
		prob[9] = 0;
		prob[10] = 0;
		prob[11] = 0;
		prob[12] = (oneMinusRSquared*(powD1 + powOneMinus2R)) / (16 * pow2);
		prob[13] = (oneMinusR*(powD1 + powOneMinus2R)*r) / (16 * pow2);
		prob[14] = ((powD1 + powOneMinus2R)*rSquared) / (16 * pow2);
		prob[15] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)) / (16 * pow2);
		prob[16] = -(oneMinusR*(-powD1 + powOneMinus2R)*r) / (16 * pow2);
		prob[17] = ((powD1 - powOneMinus2R)*rSquared) / (16 * pow2);
	}
	else
	{
		double onePlus2RInverse = 1.0 / (1 + 2 * r);
		double pow2 = std::pow(2, selfingGenerations);
		double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
		double oneMinusR = 1 - r;
		double complex1 = std::pow(1 - 2 * oneMinusR*r, selfingGenerations);
		double complex2 = -complex1 + powOneMinus2R;
		double complex3 = complex1 / pow2;
		double oneMinusRSquared = oneMinusR*oneMinusR;
		double complex4 = -2 + powOneMinus2R - 6 * r + 6 * pow2*r - 3 * powOneMinus2R*r - 4 * r * r + 2 * powOneMinus2R * r * r + complex1*(1 + r)*(1 + 2 * r);
		double onePlus2R = 1 + 2 * r;
		prob[0] = (oneMinusR*onePlus2RInverse*(-2 + complex1*onePlus2R + 2 * pow2 - powOneMinus2R - 4 * r + 2 * powOneMinus2R*r)) / (8 * pow2);
		prob[1] = -((-1 + complex1)*oneMinusR) / (24 * pow2);
		prob[2] = -((-1 + complex1)*oneMinusR) / (24 * pow2);
		prob[3] = (complex4*onePlus2RInverse) / (24 * pow2);
		prob[4] = -((-1 + complex1)*r) / (24 * pow2);
		prob[5] = (complex4*onePlus2RInverse) / (24 * pow2);
		prob[6] = -((-1 + complex1)*r) / (24 * pow2);
		prob[7] = (oneMinusRSquared*(complex1 + powOneMinus2R)) / (24 * pow2);
		prob[8] = (oneMinusR*(complex1 + powOneMinus2R)*r) / (48 * pow2);
		prob[9] = -(complex2*oneMinusRSquared) / (24 * pow2);
		prob[10] = -(complex2*oneMinusR*r) / (48 * pow2);
		prob[11] = (complex3*r*r) / 24;
		prob[12] = (oneMinusRSquared*(complex1 + powOneMinus2R)) / (24 * pow2);
		prob[13] = (oneMinusR*(complex1 + powOneMinus2R)*r) / (48 * pow2);
		prob[14] = (complex1*r*r) / (24 * pow2);
		prob[15] = -(complex2*oneMinusRSquared) / (24 * pow2);
		prob[16] = -(complex2*oneMinusR*r) / (48 * pow2);
		prob[17] = (complex1*r*r) / (24 * pow2);
	}
}
template<> void genotypeProbabilitiesWithIntercross<4, true>(std::array<double, 3>& prob, int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)/4 + (2 * r + 1 - powOneMinusR)/16)/(1 + 2 * r);
	prob[1] = prob[2] = (1 - 4*prob[0])/12;
}
template<> void genotypeProbabilitiesWithIntercross<4, false>(std::array<double, 18>& prob, int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		//There's a mistake in the mathematica derivation, so this is off by 1. I'll come back and fix it later.
		nAIGenerations = nAIGenerations - 1;
		double onePlus2RInverse = 1 / (1 + 2 * r);
		double oneMinusR = 1 - r;
		double onePlus2R = 1 + 2 * r;
		double powOneMinusR1 = std::pow(oneMinusR, nAIGenerations);
		double powOneMinusR2 = powOneMinusR1*powOneMinusR1;
		double oneMinus2R = 1 - 2 * r;
		double oneMinus2RCubed = oneMinus2R*oneMinus2R*oneMinus2R;
		double powOneMinus2R1 = std::pow(oneMinus2R, selfingGenerations);
		double quadraticPower1 = std::pow((1 + 2 * (-1 + r)*r), selfingGenerations);
		double quadratic = (3 + 4 * (-2 + r)*r), quadraticSquared = quadratic*quadratic;
		double oneMinus2RSquared = oneMinus2R*oneMinus2R;
		double oneMinus2RPow4 = oneMinus2RSquared*oneMinus2RSquared;
		double complex1 = 0.25 + powOneMinusR1*(0.75 + (-2 + r)*r);
		double complex2 = complex1 * complex1;
		double complex3 = (-1 + powOneMinusR1*(1 - 2 * r))*(-1 + powOneMinusR1*(1 - 2 * r));
		double complex4 = (-1 + oneMinus2RSquared*powOneMinusR1)*(-1 + oneMinus2RSquared*powOneMinusR1);
		double complex5 = (1 + powOneMinusR1*quadratic)*(1 + powOneMinusR1*quadratic);
		double complex6 = (0.25 - (oneMinus2RSquared*powOneMinusR1) / 4)*(0.25 - (oneMinus2RSquared*powOneMinusR1) / 4);
		double twoRMinus3 = 2 * r - 3;
		double complex7 = (1 - oneMinus2R*powOneMinusR1)*(1 - oneMinus2R*powOneMinusR1);
		double pow21 = std::pow(0.5, 3 + selfingGenerations);
		double pow22 = pow21 * 0.5;
		double pow23 = pow22 * 0.5;
		double pow24 = pow23 * 0.5;
		double pow25 = pow24 * 0.5;
		double pow26 = std::pow(2, 1 + selfingGenerations);

		prob[0] = onePlus2RInverse*pow24*(8 * onePlus2R*(-3 + pow26) + onePlus2R*quadraticPower1*(9 + 2 * powOneMinusR1*quadratic + powOneMinusR2*quadraticSquared) + 8 * powOneMinusR1*quadratic*(-1 + pow26 - oneMinus2R*powOneMinus2R1 - 2 * r));
		prob[1] = pow22*(4 + 4 * powOneMinusR1 - 8 * powOneMinusR1*r - quadraticPower1*(3 + 2 * oneMinus2R*powOneMinusR1 - oneMinus2RCubed*powOneMinusR2*twoRMinus3));
		prob[2] = pow21*(4 + 4 * oneMinus2R*oneMinusR*powOneMinusR1 + quadraticPower1*(-3 - 2 * oneMinus2R*oneMinusR*powOneMinusR1 + oneMinus2RSquared*powOneMinusR2*twoRMinus3));
		prob[3] = onePlus2RInverse*pow24*(8 * onePlus2R*(-3 + pow26) + onePlus2R*(9 - 2 * oneMinus2RSquared*powOneMinusR1 + oneMinus2RPow4*powOneMinusR2)*quadraticPower1 - 8 * oneMinus2RSquared*powOneMinusR1*(-1 + pow26 - oneMinus2R*powOneMinus2R1 - 2 * r));
		prob[4] = pow21*(4 - 4 * oneMinus2R*oneMinusR*powOneMinusR1 + (-3 + 2 * oneMinus2R*oneMinusR*powOneMinusR1 + oneMinus2RCubed*powOneMinusR2)*quadraticPower1);
		prob[5] = onePlus2RInverse*pow23*(8 * onePlus2R*(-3 + pow26) + onePlus2R*(9 - 2 * oneMinus2R*powOneMinusR1 + oneMinus2RSquared*powOneMinusR2)*quadraticPower1 - 8 * oneMinus2R*powOneMinusR1*(-1 + pow26 - oneMinus2R*powOneMinus2R1 - 2 * r));
		prob[6] = pow22*(1 - oneMinus2R*powOneMinusR1)*(4 + (-3 - oneMinus2R*powOneMinusR1)*quadraticPower1);
		prob[7] = pow25*(16 * complex6*(-powOneMinus2R1 + quadraticPower1) + complex5*(powOneMinus2R1 + quadraticPower1));
		prob[8] = pow21*(2 * oneMinus2R*oneMinusR*powOneMinus2R1*powOneMinusR1*(1 - oneMinus2R*powOneMinusR1) + (1 - oneMinus2RSquared*powOneMinusR2)*quadraticPower1);
		prob[9] = pow25*(16 * complex2*(-powOneMinus2R1 + quadraticPower1) + complex4*(powOneMinus2R1 + quadraticPower1));
		prob[10] = pow21*(-2 * oneMinus2R*oneMinusR*powOneMinus2R1*powOneMinusR1*(1 - oneMinus2R*powOneMinusR1) + (1 - oneMinus2RSquared*powOneMinusR2)*quadraticPower1);
		prob[11] = complex7*pow23*quadraticPower1;
		prob[12] = pow24*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + complex5*(powOneMinus2R1 + quadraticPower1));
		prob[13] = pow23*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + 16 * complex1*(1.0 / 4.0 - (oneMinus2RSquared*powOneMinusR1) / 4.0)*(powOneMinus2R1 + quadraticPower1));
		prob[14] = pow24*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + complex4*(powOneMinus2R1 + quadraticPower1));
		prob[15] = pow24*(16 * complex2*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
		prob[16] = pow23*(16 * complex1*(1.0 / 4.0 - (oneMinus2RSquared*powOneMinusR1) / 4)*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
		prob[17] = pow24*(16 * complex6*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
	}
	else
	{
		double oneMinusR = 1 - r;
		double powOneMinusR1 = std::pow(oneMinusR, nAIGenerations + 1);
		double powOneMinusR2 = std::pow(oneMinusR, nAIGenerations - 1);
		double powOneMinusR3 = std::pow(oneMinusR, 2 * nAIGenerations);
		double pow2 = std::pow(2, selfingGenerations);
		double quadraticPower = std::pow(1 + 2 * (r - 1) * r, selfingGenerations);
		double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
		double quadratic2 = (2 * r - 3) * (2 * r - 1);
		double quadratic1 = 3 + 4 * (r - 2) * r;
		double quadratic1Squared = quadratic1 * quadratic1;
		
		double cubicPower = -3 * oneMinusR + oneMinusR * powOneMinusR2 * quadratic1;
		cubicPower *= cubicPower;
		
		double squaredFactor1 = 1 + powOneMinusR2 * quadratic1;
		squaredFactor1 *= squaredFactor1;
		
		double oneMinusRSquared = (1 - r)*(1 - r);
		double oneMinus2R = 1 - 2 * r;
		double onePlus2R = 1 + 2 * r;
		prob[0] = (8 * oneMinusRSquared*onePlus2R*(-3 + 2 * pow2) + onePlus2R*(9 * oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2 * powOneMinusR1*quadratic2)*quadraticPower + 8 * powOneMinusR1*quadratic2*(-1 + 2 * pow2 - oneMinus2R*powOneMinus2R - 2 * r)) / (64 * oneMinusRSquared*onePlus2R*pow2);
		prob[1] = (12 * oneMinusRSquared + 4 * powOneMinusR1*quadratic2 - (9 * oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2 * powOneMinusR1*quadratic2)*quadraticPower) / (48 * oneMinusRSquared*pow2);
		prob[2] = (12*oneMinusRSquared + 4*powOneMinusR1*quadratic2 - (9*oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2*powOneMinusR1*quadratic2)*quadraticPower)/(24*oneMinusRSquared*pow2);
		prob[3] = (72*oneMinusRSquared*onePlus2R*(-3 + 2*pow2) + onePlus2R*(81*oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6*powOneMinusR1*quadratic2)*quadraticPower - 24*powOneMinusR1*quadratic2*(-1 + 2*pow2 - oneMinus2R*powOneMinus2R - 2*r))/(576*oneMinusRSquared*onePlus2R*pow2);
		prob[4] = ((-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12*oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9*r)))/(72*oneMinusRSquared*pow2);
		prob[5] = (72 * oneMinusRSquared*onePlus2R*(-3 + 2 * pow2) + onePlus2R*(81 * oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6 * powOneMinusR1*quadratic2)*quadraticPower - 24 * powOneMinusR1*quadratic2*(-1 + 2 * pow2 - oneMinus2R*powOneMinus2R - 2 * r)) / (288 * oneMinusRSquared*onePlus2R*pow2);
		prob[6] = ((-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12 * oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9 * r))) / (144 * oneMinusRSquared*pow2);
		prob[7] = (-((cubicPower*(powOneMinus2R - quadraticPower)) / oneMinusRSquared) + 9 * (powOneMinus2R + quadraticPower)*squaredFactor1) / (1152 * pow2);
		prob[8] = (-2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (72 * oneMinusRSquared*pow2);
		prob[9] = ((cubicPower*(powOneMinus2R + quadraticPower)) / oneMinusRSquared + 9 * (-powOneMinus2R + quadraticPower)*squaredFactor1) / (1152 * pow2);
		prob[10] = (2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (72 * oneMinusRSquared*pow2);
		prob[11] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
		prob[12] = (-((cubicPower*(powOneMinus2R - quadraticPower)) / oneMinusRSquared) + 9 * (powOneMinus2R + quadraticPower)*squaredFactor1) / (576 * pow2);
		prob[13] = (-2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (144 * oneMinusRSquared*pow2);
		prob[14] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
		prob[15] = ((cubicPower*(powOneMinus2R + quadraticPower)) / oneMinusRSquared + 9 * (-powOneMinus2R + quadraticPower)*squaredFactor1) / (576 * pow2);
		prob[16] = (2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (144 * oneMinusRSquared*pow2);
		prob[17] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
	}
#ifdef INTERNAL_CHECKS
	double sum = 0;
	for (int i = 0; i < 18; i++) sum += prob[i];
	if(fabs(sum - 1) > 1e-6) throw std::runtime_error("Internal error");
#endif
	prob[0] /= 4;
	prob[1] /= 16;
	prob[2] /= 32;
	prob[3] /= 4;
	prob[4] /= 32;
	prob[5] /= 8;
	prob[6] /= 16;
	prob[7] /= 4;
	prob[8] /= 32;
	prob[9] /= 4;
	prob[10] /= 32;
	prob[11] /= 8;
	prob[12] /= 8;
	prob[13] /= 16;
	prob[14] /= 8;
	prob[15] /= 8;
	prob[16] /= 16;
	prob[17] /= 8;
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<4, true>(array2<4>&data, int selfingGenerations, std::size_t nFunnels)
{
	data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25;
	data.values[0][1] = data.values[0][2] = data.values[0][3] = data.values[1][0] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[2][3] = data.values[3][0] = data.values[3][1] = data.values[3][2] = 0;
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<4, false>(array2<4>&data, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels > 1)
	{
		double pow2 = std::pow(0.5, selfingGenerations);
		data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25*(1 - pow2);
		data.values[0][1] = data.values[0][2] = data.values[0][3] = data.values[1][0] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[2][3] = data.values[3][0] = data.values[3][1] = data.values[3][2] = pow2 / 12;
	}
	else
	{
		double pow2 = std::pow(0.5, selfingGenerations);
		data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25*(1 - pow2);
		data.values[0][1] = data.values[1][0] = data.values[2][3] = data.values[3][2] = 0;
		data.values[0][2] = data.values[0][3] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[3][0] = data.values[3][1] = pow2 / 8;
	}
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<4, true>(array2<4>& data, int selfingGenerations, std::size_t nFunnels)
{
	data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25;
	data.values[0][1] = data.values[0][2] = data.values[0][3] = data.values[1][0] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[2][3] = data.values[3][0] = data.values[3][1] = data.values[3][2] = 0;
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<4, false>(array2<4>& data, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels > 1)
	{
		double pow2 = std::pow(0.5, selfingGenerations);
		data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25 - 3.0 * pow2 / 16.0;
		data.values[0][1] = data.values[0][2] = data.values[0][3] = data.values[1][0] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[2][3] = data.values[3][0] = data.values[3][1] = data.values[3][2] = pow2 / 16;
	}
	else
	{
		double pow2 = std::pow(0.5, selfingGenerations);
		data.values[0][0] = data.values[1][1] = data.values[2][2] = data.values[3][3] = 0.25 - 3.0 * pow2 / 16.0;
		data.values[0][1] = data.values[0][2] = data.values[0][3] = data.values[1][0] = data.values[1][2] = data.values[1][3] = data.values[2][0] = data.values[2][1] = data.values[2][3] = data.values[3][0] = data.values[3][1] = data.values[3][2] = pow2 / 16;
	}
}
