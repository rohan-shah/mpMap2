#include "probabilities2.h"
#include <cmath>
#include <stdexcept>
const int probabilityData<2>::intermediateProbabilitiesMask[][3] = 
		{
			{0, 1, 2},
			{1, 0, 2},
			{2, 2, 3}
		};
const int probabilityData<2>::intermediateAllelesMask[][2] = 
		{
			{0, 2},
			{2, 1}
		};
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
	//Two heterozygotes. We don't distinguish between the two different heterozygote combinations
	prob[3] = std::pow(0.5, selfingGenerations) * quadraticPower;
	//But we might want to allow us to get back to the seperate probabilities in code that calls this function, so compute the probability for state 4 by itself. 
	prob[4] = std::pow(0.5, selfingGenerations-1) * (quadraticPower - oneMinusTwoRPower);

	//The heterozygotes will be counted twice, so divide by two (or four if there are two heterozygotes). 
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
	//We don't distinguish between the two heterozygote combinations
   	prob[3] += prob[4];
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
	//Hetrozygote probabilities are multiplied by two
	data.values[1][0] = data.values[0][1] = pow2;
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
	//Hetrozygote probabilities are multiplied by two
	data.values[1][0] = data.values[0][1] = 0.5*pow2;
}
