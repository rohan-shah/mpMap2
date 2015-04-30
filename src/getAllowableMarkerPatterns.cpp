#include "getAllowableMarkerPatterns.h"
#include <stdexcept>
#include <math.h>
#include <cstring>
template<unsigned int nFounders> void getWeights(int (&weights)[4][3], int (&valuesMarker1)[nFounders], int (&valuesMarker2)[nFounders])
{
	throw std::runtime_error("Internal error");
}
template<> void getWeights<2>(int (&weights)[4][3], int (&valuesMarker1)[2], int (&valuesMarker2)[2])
{
	const int mask[2][2] =
	{
		{0, 1},
		{1, 0}
	};
	for(unsigned int founderCounter1 = 0; founderCounter1 < 2; founderCounter1++)
	{
		for(unsigned int founderCounter2 = 0; founderCounter2 < 2; founderCounter2++)
		{
			int markerCombination = (valuesMarker1[founderCounter1] << 1) + valuesMarker2[founderCounter2];
			weights[markerCombination][mask[founderCounter1][founderCounter2]] += 1;
		}
	}
}
template<> void getWeights<4>(int (&weights)[4][3], int (&valuesMarker1)[4], int (&valuesMarker2)[4])
{
	const int mask[4][4] =
	{
		{0, 1, 1, 1},
		{1, 0, 1, 1},
		{1, 1, 0, 1},
		{1, 1, 1, 0}
	};
	for(unsigned int founderCounter1 = 0; founderCounter1 < 4; founderCounter1++)
	{
		for(unsigned int founderCounter2 = 0; founderCounter2 < 4; founderCounter2++)
		{
			int markerCombination = (valuesMarker1[founderCounter1] << 1) + valuesMarker2[founderCounter2];
			weights[markerCombination][mask[founderCounter1][founderCounter2]] += 1;
		}
	}
}
template<> void getWeights<8>(int (&weights)[4][3], int (&valuesMarker1)[8], int (&valuesMarker2)[8])
{
	const int mask[8][8] =
	{
		{0, 1, 2, 2, 2, 2, 2, 2},
		{1, 0, 2, 2, 2, 2, 2, 2},
		{2, 2, 0, 1, 2, 2, 2, 2},
		{2, 2, 1, 0, 2, 2, 2, 2},
		{2, 2, 2, 2, 0, 1, 2, 2},
		{2, 2, 2, 2, 1, 0, 2, 2},
		{2, 2, 2, 2, 2, 2, 0, 1},
		{2, 2, 2, 2, 2, 2, 1, 0}
	};
	for(unsigned int founderCounter1 = 0; founderCounter1 < 8; founderCounter1++)
	{
		for(unsigned int founderCounter2 = 0; founderCounter2 < 8; founderCounter2++)
		{
			int markerCombination = (valuesMarker1[founderCounter1] << 1) + valuesMarker2[founderCounter2];
			weights[markerCombination][mask[founderCounter1][founderCounter2]] += 1;
		}
	}
}
template<unsigned int nFounders> bool isAllowable(int (&weights)[4][3])
{
	throw std::runtime_error("Internal error");
}
template<> bool isAllowable<2>(int (&weights)[4][3])
{
	return (weights[0][0] + weights[0][1] > 0) && (weights[1][0] + weights[1][1] > 0) && (weights[2][0] + weights[2][1] > 0) && (weights[3][0] + weights[3][1] > 0);
}
template<> bool isAllowable<4>(int (&weights)[4][3])
{
	return (
		weights[0][1] - 3 * weights[0][0] != 0 ||
		weights[1][1] - 3* weights[1][0] != 0 ||
		weights[2][1] - 3* weights[2][0] != 0 ||
		weights[3][1] - 3* weights[3][0] != 0
		);
}
template<> bool isAllowable<8>(int (&weights)[4][3])
{
	/*This corresponsd to a case such as:
	m1 = (0,0,1,1,0,0,1,1), m2 = (0,1,0,1,0,1,0,1)
	and also cases
	m1 = (0,0,1,1,0,0,1,1), m2 = (0,0,1,1,1,1,0,0)
*/
	bool flat = (
		((weights[0][0] == weights[0][1]) && (weights[0][2] + 3 * weights[0][1] - 9 * weights[0][0] == 0)) && 
		((weights[1][0] == weights[1][1]) && (weights[1][2] + 3 * weights[1][1] - 9 * weights[1][0] == 0)) &&
		((weights[2][0] == weights[2][1]) && (weights[2][2] + 3 * weights[2][1] - 9 * weights[2][0] == 0)) &&
		((weights[3][0] == weights[3][1]) && (weights[3][2] + 3 * weights[3][1] - 9 * weights[3][0] == 0))
		);
	if(flat) return false;
	/*The approximately uninformative case where for all of (00, 01, 10, 11), the number of recombinants is close to 7 times the number of non-recombinants. This case is also bad enough that we exclude it. This includes cases such as 
	m1 = (1,0,0,1,1,0,1,1), m2 = (0,0,1,0,0,0,1,1)
	*/
	bool almostFlat = (
		(abs(7 * weights[0][0] - weights[0][1] - weights[0][2]) <= 1) &&
		(abs(7 * weights[1][0] - weights[1][1] - weights[1][2]) <= 1) &&
		(abs(7 * weights[2][0] - weights[2][1] - weights[2][2]) <= 1) &&
		(abs(7 * weights[3][0] - weights[3][1] - weights[3][2]) <= 1)
		);
	if(almostFlat) return false;

	//if it isn't flat it could still be unidentifiable, by having a turning point for the probabilities in terms of r, in the range [0, 0.5].
	//it SEEMS as if the non-identifiability condition holds for all four probabilities (00, 30, 03 and 33) or none, so just check one. 
	if((weights[0][0] < weights[0][1] && weights[0][2] > 9*weights[0][0] - 3 * weights[0][1]) || (weights[0][0] > weights[0][1] && weights[0][2] < 9*weights[0][0] - 3 * weights[0][1]))
	{
		int numerator = 9*weights[0][0] - 3 * weights[0][1] - weights[0][2];
		int denominator = weights[0][0] - weights[0][1];
		//edge cases where the turning point is at 0 or 0.5, so there is no identifiability problem
		if(numerator == denominator || numerator == 4 * denominator) return true;
		float ratio = (float)numerator / (float)denominator;
		return ratio < 1 || ratio > 4;
	}
	return true;
}
template<unsigned int nFounders> bool isAllowableIRIP(int (&weights)[4][3])
{
	throw std::runtime_error("Internal error");
}
template<> bool isAllowableIRIP<2>(int (&weights)[4][3])
{
	return true;
}
//This is the same as for the 4-way case. 
template<> bool isAllowableIRIP<4>(int (&weights)[4][3])
{
	return (
		weights[0][1] - 3 * weights[0][0] != 0 ||
		weights[1][1] - 3* weights[1][0] != 0 ||
		weights[2][1] - 3* weights[2][0] != 0 ||
		weights[3][1] - 3* weights[3][0] != 0
		);
}
//But the 8-way case is different. Note that weights[2] and weights[1] have to be summed (there are only two different probabilities)
template<> bool isAllowableIRIP<8>(int (&weights)[4][3])
{
	/*This accounts for both flat and almost flat likelihoods. E.g, the case
	m1 = (1,0,1,0,1,0,1,0), m2 = (1,1,0,0,1,1,0,0) is flat, whereas the case
	m1 = (1,0,0,1,1,0,1,1), m2 = (0,0,1,0,0,0,1,1) is almost flat. 
	Both of these cases are caught by this one condition. */
	return !(
		(abs(7 * weights[0][0] - weights[0][1] - weights[0][2]) <= 1) &&
		(abs(7 * weights[1][0] - weights[1][1] - weights[1][2]) <= 1) &&
		(abs(7 * weights[2][0] - weights[2][1] - weights[2][2]) <= 1) &&
		(abs(7 * weights[3][0] - weights[3][1] - weights[3][2]) <= 1)
		);
}
//getAllowableMarkerPatterns passes through to a templated version because the code that needs to change is deeply nested and I don't want to duplicate the outer loop
template<unsigned int nFounders> void getAllowableMarkerPatterns(std::vector<bool>& allowableMarkerPatternsStandard, std::vector<bool>& allowableMarkerPatternsIRIP, std::map<markerEncoding, markerPatternID>& markerPatterns)
{
	std::fill(allowableMarkerPatternsStandard.begin(), allowableMarkerPatternsStandard.end(), true);
	std::fill(allowableMarkerPatternsIRIP.begin(), allowableMarkerPatternsIRIP.end(), true);
	//array to decode the marker segregation patterns into
	int valuesMarker1[nFounders], valuesMarker2[nFounders];
	//we decompose P(00), etc, into integer combinations of basis functions. This array holds those weights, first dimension is 00, 01, 10, 11, and the second dimension is the AA, AB, and AC. In the case of the 4-way AC will always be 0 because this never occurs. 
	int weights[4][3];

	for(std::map<markerEncoding, markerPatternID>::iterator markerPattern1 = markerPatterns.begin(); markerPattern1 != markerPatterns.end(); markerPattern1++)
	{
		int encodedMarkerPattern1 = markerPattern1->first;
		//We need to figure out how many alleles are present here
		int allelesPresent1[8] = {0,0,0,0,0,0,0,0};
		for(unsigned int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			int alleleValue = (encodedMarkerPattern1 & (7 << 3*founderCounter)) >> (3*founderCounter);
			allelesPresent1[alleleValue] = 1;
			valuesMarker1[founderCounter] = alleleValue;
		}
		int nAlleles1 = allelesPresent1[0] + allelesPresent1[1] + allelesPresent1[2] + allelesPresent1[3] + allelesPresent1[4] + allelesPresent1[5] + allelesPresent1[6] + allelesPresent1[7];
		//For now we don't attempt to pick up cases where multi-allelic markers lead to problems. It's acceptable NOT to touch allowableMarkerPatternsPtr in this case as we set everything to true at the start of the function. 
		if(nAlleles1 > 2) continue;
		for(std::map<markerEncoding, markerPatternID>::iterator markerPattern2 = markerPattern1; markerPattern2 != markerPatterns.end(); markerPattern2++)
		{
			memset(&(weights[0][0]), 0, sizeof(int)*12);

			int encodedMarkerPattern2 = markerPattern2->first;
			int allelesPresent2[8] = {0,0,0,0,0,0,0,0};
			for(unsigned int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				int alleleValue = (encodedMarkerPattern2 & (7 << 3*founderCounter)) >> (3*founderCounter);
				allelesPresent2[alleleValue] = 1;
				valuesMarker2[founderCounter] = alleleValue;
			}
			int nAlleles2 = allelesPresent2[0] + allelesPresent2[1] + allelesPresent2[2] + allelesPresent2[3] + allelesPresent2[4] + allelesPresent2[5] + allelesPresent2[6] + allelesPresent2[7];
			//Again, assume that multi-allelic markers are ok by default. 
			if(nAlleles2 > 2) continue;
			//Ok, now we have our pair of marker patterns, and we need to work out the form of the function that calculates the probability of a 00, 01, 11 or 10. This is going to be an integer combination of three basis functions (probabilities of AA, AB, and AC - The others such as BC, EF, etc are all just combinations of these)
			getWeights<nFounders>(weights, valuesMarker1, valuesMarker2);
			int marker1ToIndex = markerPattern1->second, marker2ToIndex = markerPattern2->second;
			bool allowableStandard;
			//non-polymorphic markers are marked as not informative
			if(nAlleles1 == 1 || nAlleles2 == 1) allowableStandard = false;
			else allowableStandard = isAllowable<nFounders>(weights); 

			allowableMarkerPatternsStandard[marker1ToIndex*markerPatterns.size() + marker2ToIndex] = allowableMarkerPatternsStandard[marker2ToIndex*markerPatterns.size() + marker1ToIndex] = allowableStandard;

			//Also check if this combination is OK for IRIP lines. 
			bool allowableIRIP;
			if(nAlleles1 == 1 || nAlleles2 == 1) allowableIRIP = false;
			else allowableIRIP = isAllowableIRIP<nFounders>(weights);
			allowableMarkerPatternsIRIP[marker1ToIndex*markerPatterns.size() + marker2ToIndex] = allowableMarkerPatternsIRIP[marker2ToIndex*markerPatterns.size() + marker1ToIndex] = allowableIRIP;
		}
	}
}
void getAllowableMarkerPatterns(std::vector<bool>& allowableMarkerPatterns, std::vector<bool>& allowableMarkerPatternsIRIP, std::map<markerEncoding, markerPatternID>& markerPatterns, unsigned int nFounders)
{
	if(nFounders == 2)
	{
		getAllowableMarkerPatterns<2>(allowableMarkerPatterns, allowableMarkerPatternsIRIP, markerPatterns);
	}
	else if(nFounders == 4)
	{
		getAllowableMarkerPatterns<4>(allowableMarkerPatterns, allowableMarkerPatternsIRIP, markerPatterns);
	}
	else if(nFounders == 8)
	{
		getAllowableMarkerPatterns<8>(allowableMarkerPatterns, allowableMarkerPatternsIRIP, markerPatterns);
	}
	else throw std::runtime_error("Internal error");
}
