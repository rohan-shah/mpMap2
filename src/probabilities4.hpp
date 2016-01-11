#ifndef PROBABILITIES_4_HEADER_GUARD
#define PROBABILITIES_4_HEADER_GUARD
#include "probabilities.hpp"
#include "probabilities2.hpp"
template<> struct probabilityData<4>
{
public:
	/*See Karl Bromans paper on intermediate generations. This mask converts allele encodings (0 - 2) into indices into 
	the array of 4 different probabilities. In terms of indices, 
	Zero = homozygote, One = other homozygote, Two = hetrozygote
	In terms of values, see Table one of the paper. Zero = first equation of table, one = second equation, etc. Note that
	we combine equations 4 and 5 into a single state. 
	*/
	static const int intermediateProbabilitiesMask[16][16];
	/*This mask takes in the two alleles at a *single* location and returns a value encoding that genotype. */
	static const int intermediateAllelesMask[4][4];
	static const int infiniteMask[4][4];
private:
	probabilityData(){}
};
const int probabilityData<4>::intermediateProbabilitiesMask[][16] = 
		{/*      0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15*/
		/*0*/	{0,  3,  10, 10, 3,  1,  6,  6,  10, 6,  2,  13, 10, 6,  13, 2},
		/*1*/	{3,  4,  9,  9,  5,  3,  7,  7,  7,  9,  13, 17, 7,  9,  17, 13},
		/*2*/	{10, 9,  12, 14, 7,  6,  14, 8,  15, 11, 10, 7,  11, 16, 9,  6},
		/*3*/	{10, 9,  14, 12, 7,  6,  8,  14, 11, 16, 6,  9,  15, 11, 7,  10},
		/*4*/	{3,  5,  7,  7,  4,  3,  9,  9,  9,  7,  13, 17, 9,  7,  17, 13},
		/*5*/	{1,  3,  6,  6,  3,  0,  10, 10, 6,  10, 2,  13, 6,  10, 13, 2},
		/*6*/	{6,  7,  14, 8,  9,  10, 12, 14, 11, 15, 10, 7,  16, 11, 9,  6},
		/*7*/	{6,  7,  8,  14, 9,  10, 14, 12, 16, 11, 6,  9,  11, 15, 7,  10},
		/*8*/	{10, 7,  15, 11, 9,  6,  11, 16, 12, 14, 10, 9,  14, 8,  7,  6},
		/*9*/	{6,  9,  11, 16, 7,  10, 15, 11, 14, 12, 10, 9,  8,  14, 7,  6},
		/*10*/	{2,  13, 10, 6,  13, 2,  10, 6,  10, 10, 0,  3,  6,  6,  3,  1},
		/*11*/	{13, 17, 7,  9,  17, 13, 7,  9,  9,  9,  3,  4,  7,  7,  5,  3},
		/*12*/	{10, 7,  11, 15, 9,  6,  16, 11, 14, 8,  6,  7,  12, 14, 9,  10},
		/*13*/	{6,  9,  16, 11, 7,  10, 11, 15, 8,  14, 6,  7,  14, 12, 9,  10},
		/*14*/	{13, 17, 9,  7,  17, 13, 9,  7,  7,  7,  3,  5,  9,  9,  4,  3},
		/*15*/	{2,  13, 6,  10, 13, 2,  6,  10, 6,  6,  1,  3,  10, 10, 3,  0}
		};
/* 0 = AA
   1 = AB
   2 = AC
   3 = AD
   4 = BA
   5 = BB
   6 = BC
   7 = BD
   8 = CA
   9 = CB
   10 = CC
   11 = CD
   12 = DA
   13 = DB
   14 = DC
   15 = DD
   */
const int probabilityData<4>::intermediateAllelesMask[][4] = 
		{
			{0,  1,  2,  3},
			{4,  5,  6,  7},
			{8,  9,  10, 11},
			{12, 13, 14, 15}
		};
const int probabilityData<4>::infiniteMask[][4] = 
		{
			{0, 1, 2, 2},
			{1, 0, 2, 2},
			{2, 2, 0, 1},
			{2, 2, 1, 0}
		};
template<> void genotypeProbabilitiesNoIntercross<4, true>(double(&prob)[nDifferentProbs], double r, int, std::size_t)
{
	prob[0] = (1-r)/(4*(1 + 2*r));
	prob[1] = prob[2] = r/(4*(1 + 2 * r));
}
template<> void genotypeProbabilitiesNoIntercross<4, false>(double(&prob)[nDifferentProbs], double r, int selfingGenerations, std::size_t nFunnels)
{
	//Some entries are divided by two. This is because they're going to be counted twice. E.g. AC, AC is going to be counted as both AC, AC and CA, CA, even though the second one is strictly speaking not allowed. 
	double twoWayProbs[nDifferentProbs];
	genotypeProbabilitiesNoIntercross<2, false>(twoWayProbs, r, selfingGenerations, nFunnels);
	//Also, some entries in the twoWayProbs have been divided, for much the same reason. But this is un-needed in this case, so reverse it.
	twoWayProbs[2] *= 2;
	twoWayProbs[3] *= 4;
	twoWayProbs[4] *= 8;
	//AA and AA
	prob[0] = 0.5*(1-r)*twoWayProbs[0];
	//AA and BB
	prob[1] = 0.5*r *twoWayProbs[0];
	//AA and CC
	prob[2] = 0.25*twoWayProbs[1];
	//AA and AB. This is not possible without intercrossing generations
	prob[3] = 0;
	//AB and AB. Not possible without intercrossing
	prob[4] = 0;
	//AB and BA. Not possible without intercrossing
	prob[5] = 0;
	//AA and BC. 
	prob[6] = 0.25 * r * twoWayProbs[2];
	//Divide this entry by two, because it's going to be counted twice.
	prob[6] /= 2;
	//AC and BA. Not possible without intercrossing
	prob[7] = 0;
	//AC and BD. Involves separating out the two hetrozygotes for the two-way (see the code that generates twoWayProbs).
	prob[8] = 0.25 * r * r * (twoWayProbs[3] - twoWayProbs[4]);
	//Divide this entry by two, because it's going to be counted twice.
	prob[8] /= 2;
	//AB and AC. Not possible without intercrossing
	prob[9] = 0;
	//AA and AC
	prob[10] = 0.25 * (1 - r) * twoWayProbs[2];
	//Divide this entry by two, because it's going to be counted twice.
	prob[10] /= 2;
	//AC and CB.  
	prob[11] = 0.25 * r * (1-r) * twoWayProbs[4];
	//Divide this entry by two, because it's going to be counted twice
	prob[11] /= 2;
	//AC and AC
	prob[12] = 0.25 * (1-r) * (1-r) * (twoWayProbs[3] - twoWayProbs[4]);
	//Divide this entry by two, because it's going to be counted twice. 
	prob[12] /= 2;
	//CC and AB. Not possible without intercrossing
	prob[13] = 0;
	//AC and BC.
	prob[14] = 0.25 * r * (1-r) * (twoWayProbs[3] - twoWayProbs[4]);
	//Divide this entry by two, because it's going to be counted twice.
	prob[14] /= 2;
	//AC and CA.
	prob[15] = 0.25 * (1 - r) * (1 - r) * twoWayProbs[4];
	//Divide this entry by two, because it's going to be counted twice.
	prob[15] /= 2;
	//AC and DB. 
	prob[16] = 0.25 * r * r * twoWayProbs[4];
	//Divide this entry by two, because it's going to be counted twice.
	prob[16] /= 2;
	//AB and CD. Not possible without intercrossing
	prob[17] = 0;
}
template<> void genotypeProbabilitiesWithIntercross<4, true>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)/4 + (2 * r + 1 - powOneMinusR)/16)/(1 + 2 * r);
	prob[1] = prob[2] = (1 - 4*prob[0])/12;
}
template<> void genotypeProbabilitiesWithIntercross<4, false>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		double onePlus2RInverse = 1/(1 + 2*r);
		double oneMinusR = 1 - r;
		double onePlus2R = 1 + 2 * r;
		double powOneMinusR1 = std::pow(oneMinusR, nAIGenerations);
		double powOneMinusR2 = powOneMinusR1*powOneMinusR1;
		double oneMinus2R = 1 - 2*r;
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
		prob[1] = onePlus2RInverse*pow24*(8 * onePlus2R*(-3 + pow26) + onePlus2R*(9 - 2 * oneMinus2RSquared*powOneMinusR1 + oneMinus2RPow4*powOneMinusR2)*quadraticPower1 - 8 * oneMinus2RSquared*powOneMinusR1*(-1 + pow26 - oneMinus2R*powOneMinus2R1 - 2 * r));
		prob[2] = onePlus2RInverse*pow23*(8 * onePlus2R*(-3 + pow26) + onePlus2R*(9 - 2 * oneMinus2R*powOneMinusR1 + oneMinus2RSquared*powOneMinusR2)*quadraticPower1 - 8 * oneMinus2R*powOneMinusR1*(-1 + pow26 - oneMinus2R*powOneMinus2R1 - 2 * r));
		prob[3] = pow22*(4 + 4 * powOneMinusR1 - 8 * powOneMinusR1*r - quadraticPower1*(3 + 2 * oneMinus2R*powOneMinusR1 - oneMinus2RCubed*powOneMinusR2*twoRMinus3));
		prob[4] = pow25*(16 * complex6*(-powOneMinus2R1 + quadraticPower1) + complex5*(powOneMinus2R1 + quadraticPower1));
		prob[5] = pow25*(16 * complex2*(-powOneMinus2R1 + quadraticPower1) + complex4*(powOneMinus2R1 + quadraticPower1));
		prob[6] = pow21*(4 - 4 * oneMinus2R*oneMinusR*powOneMinusR1 + (-3 + 2 * oneMinus2R*oneMinusR*powOneMinusR1 + oneMinus2RCubed*powOneMinusR2)*quadraticPower1);
		prob[7] = pow21*(-2 * oneMinus2R*oneMinusR*powOneMinus2R1*powOneMinusR1*(1 - oneMinus2R*powOneMinusR1) + (1 - oneMinus2RSquared*powOneMinusR2)*quadraticPower1);
		prob[8] = pow24*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + complex4*(powOneMinus2R1 + quadraticPower1));
		prob[9] = pow21*(2 * oneMinus2R*oneMinusR*powOneMinus2R1*powOneMinusR1*(1 - oneMinus2R*powOneMinusR1) + (1 - oneMinus2RSquared*powOneMinusR2)*quadraticPower1);
		prob[10] = pow21*(4 + 4 * oneMinus2R*oneMinusR*powOneMinusR1 + quadraticPower1*(-3 - 2 * oneMinus2R*oneMinusR*powOneMinusR1 + oneMinus2RSquared*powOneMinusR2*twoRMinus3));
		prob[11] = pow23*(16 * complex1*(0.25 - (oneMinus2RSquared*powOneMinusR1) / 4)*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
		prob[12] = pow24*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + complex5*(powOneMinus2R1 + quadraticPower1));
		prob[13] = pow22*(1 - oneMinus2R*powOneMinusR1)*(4 + (-3 - oneMinus2R*powOneMinusR1)*quadraticPower1);
		prob[14] = pow23*(-(complex3*(powOneMinus2R1 - quadraticPower1)) + 16 * complex1*(0.25 - (oneMinus2RSquared*powOneMinusR1) / 4)*(powOneMinus2R1 + quadraticPower1));
		prob[15] = pow24*(16 * complex2*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
		prob[16] = pow24*(16 * complex6*(-powOneMinus2R1 + quadraticPower1) + complex3*(powOneMinus2R1 + quadraticPower1));
		prob[17] = complex7*pow23*quadraticPower1;
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
		prob[1] = (72 * oneMinusRSquared*onePlus2R*(-3 + 2 * pow2) + onePlus2R*(81 * oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6 * powOneMinusR1*quadratic2)*quadraticPower - 24 * powOneMinusR1*quadratic2*(-1 + 2 * pow2 - oneMinus2R*powOneMinus2R - 2 * r)) / (576 * oneMinusRSquared*onePlus2R*pow2);
		prob[2] = (72 * oneMinusRSquared*onePlus2R*(-3 + 2 * pow2) + onePlus2R*(81 * oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6 * powOneMinusR1*quadratic2)*quadraticPower - 24 * powOneMinusR1*quadratic2*(-1 + 2 * pow2 - oneMinus2R*powOneMinus2R - 2 * r)) / (288 * oneMinusRSquared*onePlus2R*pow2);
		prob[3] = (12 * oneMinusRSquared + 4 * powOneMinusR1*quadratic2 - (9 * oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2 * powOneMinusR1*quadratic2)*quadraticPower) / (48 * oneMinusRSquared*pow2);
		prob[4] = (-((cubicPower*(powOneMinus2R - quadraticPower)) / oneMinusRSquared) + 9 * (powOneMinus2R + quadraticPower)*squaredFactor1) / (1152 * pow2);
		prob[5] = ((cubicPower*(powOneMinus2R + quadraticPower)) / oneMinusRSquared + 9 * (-powOneMinus2R + quadraticPower)*squaredFactor1) / (1152 * pow2);
		prob[6] = ((-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12 * oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9 * r))) / (72 * oneMinusRSquared*pow2);
		prob[7] = (2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (72 * oneMinusRSquared*pow2);
		prob[8] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
		prob[9] = (-2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (72 * oneMinusRSquared*pow2);
		prob[10] = (12 * oneMinusRSquared + 4 * powOneMinusR1*quadratic2 - (9 * oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2 * powOneMinusR1*quadratic2)*quadraticPower) / (24 * oneMinusRSquared*pow2);
		prob[11] = (2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (144 * oneMinusRSquared*pow2);
		prob[12] = (-((cubicPower*(powOneMinus2R - quadraticPower)) / oneMinusRSquared) + 9 * (powOneMinus2R + quadraticPower)*squaredFactor1) / (576 * pow2);
		prob[13] = ((-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12 * oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9 * r))) / (144 * oneMinusRSquared*pow2);
		prob[14] = (-2 * oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3 * oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9 * oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower) / (144 * oneMinusRSquared*pow2);
		prob[15] = ((cubicPower*(powOneMinus2R + quadraticPower)) / oneMinusRSquared + 9 * (-powOneMinus2R + quadraticPower)*squaredFactor1) / (576 * pow2);
		prob[16] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
		prob[17] = (cubicPower*quadraticPower) / (288 * oneMinusRSquared*pow2);
	}
#ifndef NDEBUG
	double sum = 0;
	for(int i = 0; i < 18; i++) sum += prob[i];
#endif
	//This is because we combined some states (see mathematica code)
	prob[0] /= 4;
	prob[1] /= 4;
	prob[2] /= 8;

	prob[3] /= 16;
	prob[4] /= 4;
	prob[5] /= 4;
	prob[6] /= 32;
	prob[7] /= 32;
	prob[8] /= 8;
	prob[9] /= 32;
	prob[10] /= 32;
	prob[11] /= 16;
	prob[12] /= 8;
	prob[13] /= 16;
	prob[14] /= 16;
	prob[15] /= 8;
	prob[16] /= 8;
	prob[17] /= 8;
}
#endif
