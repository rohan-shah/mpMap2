#include "probabilities2.h"
#include <cmath>
#include <stdexcept>
//Matrix used to turn a pair of genotypes into a probability. 
const int probabilityData<2>::intermediateProbabilitiesMask[][4] = 
		{
			{0, 2, 2, 1},
			{2, 3, 4, 2},
			{2, 4, 3, 2}, 
			{1, 2, 2, 0}
		};
//Matrix used to encode a pair of founders into a genotype.
const int probabilityData<2>::intermediateAllelesMask[][2] = 
		{
			{0, 1},
			{2, 3}
		};
//In the case of infinite generations of selfing, this table is used to turn the computed probabilities (there are two unique values) into a 2 x 2 probability matrix (containing only those two unique values)
const int probabilityData<2>::infiniteMask[][2] = 
		{
			{0, 1},
			{1, 0}
		};
template<> void genotypeProbabilitiesNoIntercross<2, true>(std::array<double, 2>& prob, double r, int, std::size_t)
{
	prob[0] = 1/(2*(1 + 2*r));
	prob[1] = r/(1 + 2 * r);
}
template<> void genotypeProbabilitiesNoIntercross<2, false>(std::array<double, 5>& prob, double r, int selfingGenerations, std::size_t)
{
	double quadraticPower = std::pow(1 - 2*r+ 2 * r * r, selfingGenerations);
	double oneMinusTwoRPower = std::pow(1 - 2 * r, selfingGenerations);
	//Probability of the same homozygote at both loci
	prob[0] = (1/(2*(1 + 2*r))) -std::pow(0.5, selfingGenerations+2)*(2 - quadraticPower + oneMinusTwoRPower*(1 - 2 *r)/ ( 1 + 2*r));
	//Different homozygotes at both
	prob[1] = (r/(1 + 2*r)) -std::pow(0.5, selfingGenerations+2)*(2 - quadraticPower - oneMinusTwoRPower*(1 - 2 *r)/ ( 1 + 2*r));
	//One homozygote, one heterozygote
	prob[2] = std::pow(0.5, selfingGenerations+1) * (1 - quadraticPower);
	//Two heterozygotes, same phase
	prob[3] = std::pow(0.5, selfingGenerations) * (quadraticPower + oneMinusTwoRPower);
	//Two heterozygotes, different phase
	prob[4] = std::pow(0.5, selfingGenerations) * (quadraticPower - oneMinusTwoRPower);

	//Computations were done with aggregated states
	prob[2] /= 2;
	prob[3] /= 4;
	prob[4] /= 4;
}
template<> void genotypeProbabilitiesWithIntercross<2, true>(std::array<double, 2>& prob, int nAIGenerations, double r, int, std::size_t)
{
	double tmp = pow(1-r, nAIGenerations - 1);
	//calculated by taking the 4-way case and setting both pairs of founders to be identical
	prob[0] = (1/(1 + 2 * r)) * ((1-r)*tmp/2 + (2*r + 1 - tmp) /4);
	prob[1] = (1 - prob[0]*2)/2;
}
template<> void genotypeProbabilitiesWithIntercross<2, false>(std::array<double, 5>& prob, int nAIGenerations, double r, int selfingGenerations, std::size_t)
{
	double pow2 = std::pow(2, selfingGenerations);
	double powOneMinusR = std::pow(1 - r, 1 + nAIGenerations);
	double powOneMinusRSquared = std::pow(1 - r, 2*nAIGenerations);
	double quadraticPower = std::pow(1 + 2* (-1 + r)* r, selfingGenerations);
	double powOneMinus2R = std::pow(1 - 2*r, 1 + selfingGenerations);
	double oneMinusTwoRSquared = (1 - 2*r)*(1 - 2*r);
	double onePlus2R = 1 + 2*r;
	double oneMinusRSquared = (1 - r) * (1 - r);
	//Same homozygote at both loci
	prob[0] = 
	0.5*(4*oneMinusRSquared*onePlus2R*(-1 + pow2) 
		+ onePlus2R*(oneMinusRSquared + oneMinusTwoRSquared*powOneMinusRSquared)*quadraticPower 
		- 2*(2*pow2 - powOneMinus2R)*powOneMinusR*(-1 + 2*r)
	)/(8*onePlus2R*pow2*oneMinusRSquared);

	prob[1] = 
  	0.5*(4*oneMinusRSquared*onePlus2R*(-1 + pow2) 
  		+ onePlus2R*(oneMinusRSquared + oneMinusTwoRSquared*powOneMinusRSquared)*quadraticPower 
  		+ 2*(2*pow2 - powOneMinus2R)*powOneMinusR*(-1 + 2*r)
  	)/(8*onePlus2R*pow2*oneMinusRSquared);

 	prob[2] = 
 	0.125*(2*oneMinusRSquared 
 		- (oneMinusRSquared + oneMinusTwoRSquared*powOneMinusRSquared)*quadraticPower
	)/(2*pow2*oneMinusRSquared);

	prob[3] = 
	0.25*(2*powOneMinus2R*powOneMinusR 
		+ (oneMinusRSquared + oneMinusTwoRSquared*powOneMinusRSquared)*quadraticPower
	)/(8*pow2*oneMinusRSquared);

   	prob[4] = 
	0.25*(-2*powOneMinus2R*powOneMinusR 
   		+ (oneMinusRSquared + oneMinusTwoRSquared*powOneMinusRSquared)*quadraticPower
   	)/(8*pow2*oneMinusRSquared);
	prob[3] *= 2;
	prob[4] *= 2;
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<2, true>(array2<2>&data, int selfingGenerations, std::size_t nFunnels)
{
	data.values[0][0] = data.values[1][1] = 0.5;
	data.values[1][0] = data.values[0][1] = 0;
}
template<> void singleLocusGenotypeProbabilitiesNoIntercross<2, false>(array2<2>&data, int selfingGenerations, std::size_t nFunnels)
{
	double pow2 = std::pow(0.5, selfingGenerations);
	data.values[0][0] = data.values[1][1] = 0.5 * (1 - pow2);
	data.values[1][0] = data.values[0][1] = 0.5 * pow2;
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<2, true>(array2<2>& data, int selfingGenerations, std::size_t nFunnels)
{
	data.values[0][0] = data.values[1][1] = 0.5;
	data.values[1][0] = data.values[0][1] = 0;
}
template<> void singleLocusGenotypeProbabilitiesWithIntercross<2, false>(array2<2>& data, int selfingGenerations, std::size_t nFunnels)
{
	double pow2 = std::pow(0.5, selfingGenerations);
	data.values[0][0] = data.values[1][1] = 0.25 * (2 - pow2);
	data.values[1][0] = data.values[0][1] = 0.25*pow2;
}
