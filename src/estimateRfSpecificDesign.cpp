#include "estimateRfSpecificDesign.h"
#include <math.h>
#include "intercrossingAndSelfingGenerations.h"
#include "estimateRfCheckFunnels.h"
#include <map>
#include <set>
#include "recodeFoundersFinalsHets.h"
#include "markerPatternsToUniqueValues.h"
#include "funnelsToUniqueValues.h"
#include <memory>
struct rfhaps_internal_args
{
	rfhaps_internal_args(std::vector<double>& lineWeights, std::vector<double>& recombinationFractions)
	: recombinationFractions(recombinationFractions), lineWeights(lineWeights)
	{}
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	Rcpp::S4 hetData;
	std::vector<double>& recombinationFractions;
	
	std::vector<int> intercrossingGenerations;
	//A vector where entry i contains the markerPatternID identifying the segregation pattern of marker number i. Contains one entry per marker. 
	std::vector<markerPatternID> markerPatternIDs;
	std::vector<double>& lineWeights;
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i
	std::vector<markerEncoding> markerEncodings;
	bool hasAI;
	int marker1Start, marker1End;
	int marker2Start, marker2End;
	//maximum number of marker alleles present
	int maxAlleles;
	//zeroed by the calling code. Must be added to, not overwritten. 
	double* result;
	std::vector<funnelID> funnelIDs;
	std::vector<funnelEncoding> funnelEncodings;
};
template<int n> struct arrayType
{
public:
	arrayType()
	{}
	double values[n][n];
};
//Note that if we change the mask[8][8] values of 2 to 1 we get mask4 in the first 4x4 block. 
//const int mask4[4][4] = {{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
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
//Templated function to work out the two-point probabilities with the given recombination fraction (and number of AI generations). Templating allows the number of founders to be a compile-time constant
template<int nFounders> void genotypeProbabilitiesNoIntercross(double (&prob)[3], double recombinationFraction);
template<int nFounders> void genotypeProbabilitiesWithIntercross(double (&prob)[3], int nAIGenarations, double recombinationFraction);
//There are only really two probability values for the 4-way design, but if we put in three we can use the same mask as for the 8-way case
template<> void genotypeProbabilitiesNoIntercross<4>(double (&prob)[3], double r)
{
	prob[0] = (1-r)/(4+8*r);
	prob[1] = prob[2] = r/(4+8*r);
}
template<> void genotypeProbabilitiesNoIntercross<8>(double (&prob)[3], double r)
{
	prob[0] = (1-r)*(1-r)/(8+16*r);
	prob[1] = r*(1-r)/(8+16*r);
	prob[2] = r/(16+32*r);
}
template<> void genotypeProbabilitiesNoIntercross<2>(double (&prob)[3], double r)
{
	prob[0] = 1/(2*(1 + 2*r));
	prob[1] = r/(1 + 2 * r);
}
template<> void genotypeProbabilitiesWithIntercross<2>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations - 1);
	//calculated by taking the 4-way case and setting both pairs of founders to be identical
	prob[0] = (1/(1 + 2 * r)) * ((1-r)*tmp/2 + (2*r + 1 - tmp) /4);
	prob[1] = (1 - prob[0]*2)/2;
}
template<> void genotypeProbabilitiesWithIntercross<4>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	//prob[0] = (pow(1-r, 1+nAIGenerations)/4+(2*r+1-pow(1-r, nAIGenerations-1))/16)/(1+2*r); 
	prob[0] = (tmp *(1-r)*(1-r)/4 + (2*r + 1 - tmp)/16)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 4 * prob[0]) / 12;
}
template<> void genotypeProbabilitiesWithIntercross<8>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	prob[0] = (tmp *(1-r)*(1-r)*(1-r)/8 + (2*r + 1 - tmp)/64)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 8 * prob[0]) / 56;
}
template<int nFounders> void genotypeProbabilitiesNoIntercross(arrayType<nFounders>& expandedProbabilities, double r)
{
	double probabilities[3];
	genotypeProbabilitiesNoIntercross<nFounders>(probabilities, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities.values[i][j] = probabilities[mask[i][j]];
		}
	}
}
template<int nFounders> void genotypeProbabilitiesWithIntercross(double (&expandedProbabilities)[nFounders][nFounders], int nAIGenerations, double r)
{
	double probabilities[3];
	genotypeProbabilitiesWithIntercross<nFounders>(probabilities, nAIGenerations, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities[i][j] = probabilities[mask[i][j]];
		}
	}
}
template<int maxAlleles> struct perRecombinationFractionData
{
	std::vector<arrayType<maxAlleles> > perFunnelData;
	std::vector<arrayType<maxAlleles> > perAIGenerationData;
};
template<int maxAlleles> struct singleMarkerPairData
{
public:
	singleMarkerPairData(int nRecombLevels)
		:allowableFunnel(true), allowableAI(true)
	{
		perRecombData = new perRecombinationFractionData<maxAlleles>[nRecombLevels];
	}
	singleMarkerPairData(const singleMarkerPairData& other)
	: perRecombData(other.perRecombData), allowableFunnel(other.allowableFunnel), allowableAI(other.allowableAI)
	{}
	perRecombinationFractionData<maxAlleles>* perRecombData;
	bool allowableFunnel;
	bool allowableAI;
};
typedef std::pair<markerPatternID, markerPatternID> markerPair;
template<int maxAlleles> class allMarkerPairData : public std::map<markerPair, singleMarkerPairData<maxAlleles> >
{
public:
	typedef typename std::map<markerPair, singleMarkerPairData<maxAlleles> > parent;
	~allMarkerPairData()
	{
		for(typename parent::iterator i = parent::begin(); i != parent::end(); i++)
		{
			delete[] i->second.perRecombData;
		}
	}
};
template<int maxAlleles> struct constructLookupTableArgs
{
public:
	constructLookupTableArgs(allMarkerPairData<maxAlleles>& computedContributions)
		: computedContributions(computedContributions)
	{}
	allMarkerPairData<maxAlleles>& computedContributions;
	std::vector<markerEncoding>* markerEncodings;
	std::vector<funnelEncoding>* funnelEncodings;
	std::vector<double>* recombinationFractions;
	std::vector<int>* intercrossingGenerations;
};
template<int nFounders, int maxAlleles, bool infiniteSelfing> void constructLookupTable(constructLookupTableArgs<maxAlleles>& args)
{
	int nMarkerPatternIDs = args.markerEncodings->size();
	int nRecombLevels = args.recombinationFractions->size();
	int nDifferentFunnels = args.funnelEncodings->size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations->begin(), args.intercrossingGenerations->end());

	//Only compute the haplotype probabilities once
	std::vector<arrayType<nFounders> > funnelHaplotypeProbabilities(nRecombLevels);
	for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
	{
		genotypeProbabilitiesNoIntercross<nFounders>(funnelHaplotypeProbabilities[recombCounter], (*args.recombinationFractions)[recombCounter]);
	}
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static, 1)
#endif
	//This is a big chunk of code, but does NOT grow with problem size (number of markers, number of lines). Well, it grows but to some fixed limit. 
	for(int firstPattern = 0; firstPattern < nMarkerPatternIDs; firstPattern++)
	{
		int firstMarkerEncoding = (*args.markerEncodings)[firstPattern];
		int firstMarkerPattern[nFounders];
		for(int i = 0; i < nFounders; i++)
		{
			firstMarkerPattern[i] = ((firstMarkerEncoding & (7 << (3*i))) >> (3*i));
		}
		//marker alleles have been encoded so they're of the form [0, nAlleles), so can just look for max value
		int nFirstMarkerAlleles = *std::max_element(firstMarkerPattern, firstMarkerPattern+nFounders);
		for(int secondPattern = 0; secondPattern < nMarkerPatternIDs; secondPattern++)
		{
			int secondMarkerEncoding = (*args.markerEncodings)[secondPattern];
			int secondMarkerPattern[nFounders];
			for(int i = 0; i < nFounders; i++)
			{
				secondMarkerPattern[i] = ((secondMarkerEncoding & (7 << (3*i))) >> (3*i));
			}
			int nSecondMarkerAlleles = *std::max_element(secondMarkerPattern, secondMarkerPattern+nFounders);			
			
			markerPair currentPair(firstPattern, secondPattern);
			singleMarkerPairData<maxAlleles> thisMarkerPairData(nRecombLevels);
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				arrayType<nFounders>& funnelHaplotypeProbabilitiesThisRecomb = funnelHaplotypeProbabilities[recombCounter];

				double recombFraction = (*args.recombinationFractions)[recombCounter];
				perRecombinationFractionData<maxAlleles>& currentRecomb = thisMarkerPairData.perRecombData[recombCounter];
				currentRecomb.perFunnelData.resize(nDifferentFunnels);
				currentRecomb.perAIGenerationData.resize(maxAIGenerations);
				for(int intercrossingGeneration = 1; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
				{
					arrayType<maxAlleles>& markerProbabilities = currentRecomb.perAIGenerationData[intercrossingGeneration-1];
					memset(&markerProbabilities, 0, sizeof(arrayType<maxAlleles>));
					double haplotypeProbabilities[nFounders][nFounders];
					genotypeProbabilitiesWithIntercross<nFounders>(haplotypeProbabilities, intercrossingGeneration, recombFraction);
					for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
					{
						for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
						{
							if(firstMarkerPattern[firstFounder] == firstMarkerValue)
							{
								for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
								{
									for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
									{
										 if(secondMarkerPattern[secondFounder] == secondMarkerValue)
										 {
											markerProbabilities.values[firstMarkerValue][secondMarkerValue] += haplotypeProbabilities[firstFounder][secondFounder];
										 }
									 }
								}
							}
						}
					}
					//now take logs of every value in markerProbabilities
					for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
					{
						for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
						{
							markerProbabilities.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilities.values[firstMarkerValue][secondMarkerValue]);
						}
					}
				}
				for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
				{
					arrayType<maxAlleles>& markerProbabilitiesThisFunnel = currentRecomb.perFunnelData[funnelCounter];
					memset(&markerProbabilitiesThisFunnel, 0, sizeof(arrayType<maxAlleles>));
					funnelEncoding enc = (*args.funnelEncodings)[funnelCounter];
					int funnel[8];
					for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
					{
						funnel[founderCounter] = ((enc & (7 << (3*founderCounter))) >> (3*founderCounter));
					}
					for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
					{
						//firstFounder is the index of the founder within the current funnel
						for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
						{
							if(firstMarkerPattern[funnel[firstFounder]] == firstMarkerValue)
							{
								for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
								{
									for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
									{
										 if(secondMarkerPattern[funnel[secondFounder]] == secondMarkerValue)
										 {
											markerProbabilitiesThisFunnel.values[firstMarkerValue][secondMarkerValue] += funnelHaplotypeProbabilitiesThisRecomb.values[firstFounder][secondFounder];
										 }
									 }
								}
							}
						}
					}
					//now take logs of every value in markerProbabilitiesThisFunnel
					for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
					{
						for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
						{
							markerProbabilitiesThisFunnel.values[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisFunnel.values[firstMarkerValue][secondMarkerValue]);
						}
					}
				}
			}
			//Now we see if these markers are informative with either 0 AIC generations or some non-zero number
#ifdef USE_OPENMP
			#pragma omp critical
#endif
			{
				args.computedContributions.insert(make_pair(currentPair, thisMarkerPairData));
			}
		}
	}
}
template<int nFounders, int maxAlleles, bool infiniteSelfing> bool estimateRfSpecificDesign(rfhaps_internal_args& args)
{
	int nFinals = args.finals.nrow(), nRecombLevels = args.recombinationFractions.size();
	int nDifferentFunnels = args.funnelEncodings.size();
	int marker2RangeSize = args.marker2End - args.marker2Start;
	int maxStart = std::max(args.marker1Start, args.marker2Start), minEnd = std::min(args.marker1End, args.marker2End);
	std::vector<double>& lineWeights = args.lineWeights;
	Rcpp::List finalDimNames = args.finals.attr("dimnames");
	Rcpp::CharacterVector finalNames = finalDimNames[0];

	int nMarkerPatternIDs = args.markerEncodings.size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations.begin(), args.intercrossingGenerations.end());

	//This is basically just a huge lookup table
	allMarkerPairData<maxAlleles> computedContributions;
	
	constructLookupTableArgs<maxAlleles> lookupArgs(computedContributions);
	lookupArgs.markerEncodings = &args.markerEncodings;
	lookupArgs.recombinationFractions = &args.recombinationFractions;
	lookupArgs.funnelEncodings = &args.funnelEncodings;
	lookupArgs.intercrossingGenerations = &args.intercrossingGenerations;
	constructLookupTable<nFounders, maxAlleles, infiniteSelfing>(lookupArgs);

#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static, 1)
#endif
	//This set of loops DOES grow with problem size. 
	for(int markerCounter1 = args.marker1Start; markerCounter1 < args.marker1End; markerCounter1++)
	{
		int markerPatternID1 = args.markerPatternIDs[markerCounter1];
		for(int markerCounter2 = args.marker2Start; markerCounter2 < args.marker2End; markerCounter2++)
		{
			//For some bits in the lower triangle we have already calculated a corresponding bit in the upper triangle, so don't recalculate these. This means we have some values in args.result which are never set, but later code is aware of this
			if(markerCounter2 >= maxStart && markerCounter2 < minEnd && markerCounter1 >= maxStart && markerCounter1 < minEnd && markerCounter2 < markerCounter1) continue;
			int markerPatternID2 = args.markerPatternIDs[markerCounter2];

			typename allMarkerPairData<maxAlleles>::iterator patternData = computedContributions.find(markerPair(markerPatternID1, markerPatternID2));
			if(patternData == computedContributions.end()) throw std::runtime_error("Internal error");
			singleMarkerPairData<maxAlleles> markerPairData = patternData->second;
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				perRecombinationFractionData<maxAlleles>& perRecombLevelData = markerPairData.perRecombData[recombCounter];
				for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
				{
					int marker1Value = args.finals(finalCounter, markerCounter1);
					int marker2Value = args.finals(finalCounter, markerCounter2);
					if(marker1Value != NA_INTEGER && marker2Value != NA_INTEGER)
					{
						double contribution = 0;
						bool allowable = false;
						int intercrossingGenerations = args.intercrossingGenerations[finalCounter];
						if(intercrossingGenerations == 0 && markerPairData.allowableFunnel)
						{
							funnelID currentLineFunnelID = args.funnelIDs[finalCounter];
							arrayType<maxAlleles>& perMarkerGenotypeValues = perRecombLevelData.perFunnelData[currentLineFunnelID];
							contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
						}
						else if(intercrossingGenerations > 0 && markerPairData.allowableAI)
						{
							arrayType<maxAlleles>& perMarkerGenotypeValues = perRecombLevelData.perAIGenerationData[intercrossingGenerations-1];
							contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
						}
						//We get an NA from trying to take the logarithm of zero - That is, this parameter is completely impossible for the given data, so put in -Inf
						if(contribution != contribution) args.result[(long)(markerCounter1 - args.marker1Start) *(long)nRecombLevels*(long)marker2RangeSize + (long)(markerCounter2-args.marker2Start) * (long)nRecombLevels + (long)recombCounter] = -std::numeric_limits<double>::infinity();
						else if(contribution != 0 && allowable) args.result[(long)(markerCounter1 - args.marker1Start) *(long)nRecombLevels*(long)marker2RangeSize + (long)(markerCounter2-args.marker2Start) * (long)nRecombLevels + (long)recombCounter] += lineWeights[finalCounter] * contribution;
					}
				}
			}
		}
	}
	return true;
}
template<int nFounders, int maxAlleles> bool estimateRfSpecificDesignInternal2(rfhaps_internal_args& args)
{
		bool infiniteSelfing = Rcpp::as<std::string>(args.pedigree.slot("selfing")) == "infinite";
		if(infiniteSelfing)
		{
			return estimateRfSpecificDesign<nFounders, maxAlleles, true>(args);
		}
		else return estimateRfSpecificDesign<nFounders, maxAlleles, false>(args);
}
//here we transfer maxAlleles over to the templated parameter section - This can make a BIG difference to memory usage if this is smaller, and it's going into a type so it has to be templated.
template<int nFounders> bool estimateRfSpecificDesignInternal1(rfhaps_internal_args& args)
{
	switch(args.maxAlleles)
	{
		case 1:
			return estimateRfSpecificDesignInternal2<nFounders, 1>(args);
		case 2:
			return estimateRfSpecificDesignInternal2<nFounders, 2>(args);
		case 3:
			return estimateRfSpecificDesignInternal2<nFounders, 3>(args);
		case 4:
			return estimateRfSpecificDesignInternal2<nFounders, 4>(args);
		case 5:
			return estimateRfSpecificDesignInternal2<nFounders, 5>(args);
		case 6:
			return estimateRfSpecificDesignInternal2<nFounders, 6>(args);
		case 7:
			return estimateRfSpecificDesignInternal2<nFounders, 7>(args);
		case 8:
			return estimateRfSpecificDesignInternal2<nFounders, 8>(args);
		default:
			throw std::runtime_error("Internal error");
	}
}
bool estimateRfSpecificDesign(estimateRfSpecificDesignArgs& args)
{
	//work out the number of intercrossing generations
	int nFounders = args.founders.nrow(), nFinals = args.finals.nrow(), nMarkers = args.finals.ncol();
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(args.pedigree, args.finals, nFounders, intercrossingGenerations, selfingGenerations);
	bool hasAIC = *std::max(intercrossingGenerations.begin(), intercrossingGenerations.end()) > 0;

	//Check that everything has a unique funnel - For the case of the lines which are just selfing, we just check that one funnel. For AIC lines, we check the funnels of all the parent lines
	std::vector<funnelType> allFunnels;
	{
		std::vector<std::string> funnelWarnings, funnelErrors;
		estimateRfCheckFunnels(args.finals, args.founders, args.pedigree, intercrossingGenerations, funnelWarnings, funnelErrors, allFunnels);
		std::size_t warningIndex = 0;
		for(;warningIndex < funnelWarnings.size() && warningIndex < 6;warningIndex++)
		{
			Rprintf(funnelWarnings[warningIndex].c_str());
		}
		if(funnelWarnings.size() > 6)
		{
			Rprintf("Supressing further funnel warnings");
		}
	}
	//re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles
	//We do this to make it easier to identify markers with identical segregation patterns. recodedFounders = column major matrix
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	Rcpp::List recodedHetData(nMarkers);
	recodedHetData.attr("names") = args.hetData.attr("names");
	recodedFinals.attr("dimnames") = args.finals.attr("dimnames");
	
	recodeDataStruct recoded;
	recoded.recodedFounders = recodedFounders;
	recoded.recodedFinals = recodedFinals;
	recoded.founders = args.founders;
	recoded.finals = args.finals;
	recoded.hetData = args.hetData;
	recoded.recodedHetData = recodedHetData;
	recodeFoundersFinalsHets(recoded);
	unsigned int maxAlleles = recoded.maxAlleles;
	if(maxAlleles > 8) 
	{
		args.error = "Internal error - Cannot have more than eight alleles per marker";
		return false;
	}
	//map containing encodings of all the marker segregation patterns, and an associated marker ID (we can't use the encoding as an index because they'll jump around a lot and could be quite big values sometimes I think). Marker IDs are guaranteed to be contiguous numbers starting from 0 - So the set of all valid [0, markerPatterns.size()]. 
	//Note that markerEncoding and markerPatternID are defined in unitTypes.hpp. They're just integers (and automatically convertible to integers), but they're represented by different types - This stops us from confusing or accidentally interchanging them.
	std::map<markerEncoding, markerPatternID> markerPatterns;
	//A vector where entry i contains the markerPatternID identifying the segregation pattern of marker number i. Contains one entry per marker - but the same value will appear multiple times because many markers will have the same segregation pattern - which is the whole point of this set-up
	std::vector<markerPatternID> markerPatternIDs;
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i. Note that this vector has a *different* number of entries to the markerPatternIDs vector.
	std::vector<markerEncoding> markerEncodings;
	markerPatternsToUniqueValues(markerPatterns, markerPatternIDs, markerEncodings, nFounders, nMarkers, recodedFounders);
	
	//map containing encodings of the funnels involved in the experiment (as key), and an associated unique index (again, using the encoded values directly is no good because they'll be all over the place). Unique indices are contiguous again.
	std::map<funnelEncoding, funnelID> funnelTranslation;
	//vector giving the funnel ID for each individual
	std::vector<funnelID> funnelIDs;
	//vector giving the encoded value for each individual
	std::vector<funnelEncoding> funnelEncodings;
	funnelsToUniqueValues(funnelTranslation, funnelIDs, funnelEncodings, allFunnels, nFounders);
	
	rfhaps_internal_args internal_args(args.lineWeights, args.recombinationFractions);
	internal_args.founders = recodedFounders;
	internal_args.finals = recodedFinals;
	internal_args.pedigree = args.pedigree;
	internal_args.hetData = recodedHetData;
	internal_args.intercrossingGenerations.swap(intercrossingGenerations);
	internal_args.markerPatternIDs.swap(markerPatternIDs);
	internal_args.markerEncodings.swap(markerEncodings);
	internal_args.hasAI = hasAIC;
	internal_args.marker1Start = args.marker1Start;
	internal_args.marker1End = args.marker1End;
	internal_args.marker2Start = args.marker2Start;
	internal_args.marker2End = args.marker2End;
	internal_args.maxAlleles = maxAlleles;
	internal_args.result = args.result;
	internal_args.funnelIDs.swap(funnelIDs);
	internal_args.funnelEncodings.swap(funnelEncodings);
	if(nFounders == 2)
	{
		return estimateRfSpecificDesignInternal1<2>(internal_args);
	}
	else if(nFounders == 4)
	{
		return estimateRfSpecificDesignInternal1<4>(internal_args);
	}
	else if(nFounders == 8)
	{
		return estimateRfSpecificDesignInternal1<8>(internal_args);
	}
	else
	{
		Rprintf("Number of founders must be 2, 4 or 8\n");
		return false;
	}
	return true;
}
