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
		/*4*/	{3,  5,  7,  7,  4,  3,  9,  9,  9,  7,  13, 17, 9,  9,  17, 13},
		/*5*/	{1,  3,  6,  6,  3,  0,  10, 10, 6,  10, 2,  13, 6,  10, 13, 2},
		/*6*/	{6,  7,  14, 8,  9,  10, 12, 14, 11, 15, 10, 7,  16, 11, 9,  6},
		/*7*/	{6,  7,  8,  14, 9,  10, 14, 12, 16, 11, 6,  9,  11, 15, 7,  10},
		/*8*/	{10, 7,  15, 11, 9,  6,  11, 16, 12, 14, 10, 9,  14, 8,  7,  6},
		/*9*/	{6,  9,  11, 16, 7,  10, 15, 11, 14, 12, 10, 9,  8,  14, 7,  6},
		/*10*/	{2,  13, 10, 6,  13, 2,  10, 6,  10, 10, 0,  3,  6,  6,  3,  1},
		/*11*/	{13, 17, 7,  9,  17, 13, 7,  9,  9,  9,  3,  4,  7,  7,  5,  3},
		/*12*/	{10, 7,  11, 15, 9,  6,  16, 11, 14, 8,  6,  7,  12, 14, 9,  10},
		/*13*/	{6,  9,  16, 11, 9,  10, 11, 15, 8,  14, 6,  7,  14, 12, 9,  10},
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
template<> void genotypeProbabilitiesNoIntercross<4, true>(double (&prob)[nDifferentProbs], double r, int)
{
	prob[0] = (1-r)/(4*(1 + 2*r));
	prob[1] = prob[2] = r/(4*(1 + 2 * r));
}
template<> void genotypeProbabilitiesNoIntercross<4, false>(double (&prob)[nDifferentProbs], double r, int selfingGenerations)
{
	//Some entries are divided by two. This is because they're going to be counted twice. E.g. AC, AC is going to be counted as both AC, AC and CA, CA, even though the second one is strictly speaking not allowed. 
	double twoWayProbs[nDifferentProbs];
	genotypeProbabilitiesNoIntercross<2, false>(twoWayProbs, r, selfingGenerations);
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
	//AC and CB. Not possible without intercrossing (But possible further down with BC instead of CB). 
	prob[11] = 0;
	//AC and AC
	prob[12] = 0.25 * (1-r) * (1-r) * (twoWayProbs[3] - twoWayProbs[4]);
	//Divide this entry by two, because it's going to be counted twice. 
	prob[12] /= 2;
	//CC and AB. Not possible without intercrossing
	prob[13] = 0;
	//AC and BC.
	prob[14] = 0.25 * r * (1-r) * (twoWayProbs[3] - twoWayProbs[4]);
	//Divide this entry by two, because it's going to be counted twice
	prob[14] /= 2;
	//AC and CA. Not possible without intercrossing (But possible with AC instaed of CA). 
	prob[15] = 0;
	//AC and DB. Not possible without intercrossing (But possible with BD instead of DB). 
	prob[16] = 0;
	//AB and CD. Not possible without intercrossing
	prob[17] = 0;
}
template<> void genotypeProbabilitiesWithIntercross<4, true>(double (&prob)[nDifferentProbs], int nAIGenerations, double r, int)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)/4 + (2 * r + 1 - powOneMinusR)/16)/(1 + 2 * r);
	prob[1] = prob[2] = (1 - 4*prob[0])/12;
}
template<> void genotypeProbabilitiesWithIntercross<4, false>(double (&prob)[nDifferentProbs], int nAIGenerations, double r, int selfingGenerations)
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
	prob[0] = (8*oneMinusRSquared*onePlus2R*(-3 + 2*pow2) + onePlus2R*(9*oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2*powOneMinusR1*quadratic2)*quadraticPower + 8*powOneMinusR1*quadratic2*(-1 + 2*pow2 - oneMinus2R*powOneMinus2R - 2*r))/(64*onePlus2R*pow2); 
	prob[1] = (72*oneMinusRSquared*onePlus2R*(-3 + 2*pow2) + onePlus2R*(81*oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6*powOneMinusR1*quadratic2)*quadraticPower - 24*powOneMinusR1*quadratic2*(-1 + 2*pow2 - oneMinus2R*powOneMinus2R - 2*r))/(576*onePlus2R*pow2);
	prob[2] = (72*oneMinusRSquared*onePlus2R*(-3 + 2*pow2) + onePlus2R*(81*oneMinusRSquared + powOneMinusR3*quadratic1Squared - 6*powOneMinusR1*quadratic2)*quadraticPower - 24*powOneMinusR1*quadratic2*(-1 + 2*pow2 - oneMinus2R*powOneMinus2R - 2*r))/(288*onePlus2R*pow2);
	prob[3] = (12*oneMinusRSquared + 4*powOneMinusR1*quadratic2 - (9*oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2*powOneMinusR1*quadratic2)*quadraticPower)/(48*pow2);
	prob[4] = (oneMinusRSquared*(-((cubicPower*(powOneMinus2R - quadraticPower))/oneMinusRSquared) + 9*(powOneMinus2R + quadraticPower)*squaredFactor1))/(1152*pow2);
	prob[5] = (oneMinusRSquared*((cubicPower*(powOneMinus2R + quadraticPower))/oneMinusRSquared + 9*(-powOneMinus2R + quadraticPower)*squaredFactor1))/(1152*pow2);
	prob[6] = ((-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12*oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9*r)))/(72*pow2);
	prob[7] = (2*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9*oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower)/(72*pow2);
	prob[8] = (cubicPower*quadraticPower)/(288*pow2);
	prob[9] = (-2*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9*oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower)/(72*pow2);
	prob[10] = (12*oneMinusRSquared + 4*powOneMinusR1*quadratic2 - (9*oneMinusRSquared + powOneMinusR3*quadratic1Squared + 2*powOneMinusR1*quadratic2)*quadraticPower)/(24*pow2);
	prob[11] = (2*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9*oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower)/(144*pow2);
	prob[12] = (oneMinusRSquared*(-((cubicPower*(powOneMinus2R - quadraticPower))/oneMinusRSquared) + 9*(powOneMinus2R + quadraticPower)*squaredFactor1))/(576*pow2);
	prob[13] = ((-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1)*(-12*oneMinusR + quadraticPower*(9 + oneMinusR*powOneMinusR2*quadratic1 - 9*r)))/(144*pow2);
	prob[14] = (-2*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1*(-3*oneMinusR + oneMinusR*powOneMinusR2*quadratic1) + (9*oneMinusRSquared - powOneMinusR3*quadratic1Squared)*quadraticPower)/(144*pow2);
	prob[15] = (oneMinusRSquared*((cubicPower*(powOneMinus2R + quadraticPower))/oneMinusRSquared + 9*(-powOneMinus2R + quadraticPower)*squaredFactor1))/(576*pow2);
	prob[16] = (cubicPower*quadraticPower)/(288*pow2);
	prob[17] = (cubicPower*quadraticPower)/(288*pow2);
}
#endif
