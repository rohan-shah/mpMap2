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
#include "constructLookupTable.hpp"
struct rfhaps_internal_args
{
	rfhaps_internal_args(std::vector<double>& lineWeights, std::vector<double>& recombinationFractions)
	: recombinationFractions(recombinationFractions), lineWeights(lineWeights)
	{}
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	std::vector<double>& recombinationFractions;
	
	std::vector<int> intercrossingGenerations;
	std::vector<int> selfingGenerations;
	std::vector<double>& lineWeights;
	markerPatternsToUniqueValuesArgs markerPatternData;
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
template<int nFounders, int maxAlleles, bool infiniteSelfing> bool estimateRfSpecificDesign(rfhaps_internal_args& args)
{
	int nFinals = args.finals.nrow(), nRecombLevels = args.recombinationFractions.size();
	int nDifferentFunnels = args.funnelEncodings.size();
	int marker2RangeSize = args.marker2End - args.marker2Start;
	int maxStart = std::max(args.marker1Start, args.marker2Start), minEnd = std::min(args.marker1End, args.marker2End);
	std::vector<double>& lineWeights = args.lineWeights;
	Rcpp::List finalDimNames = args.finals.attr("dimnames");
	Rcpp::CharacterVector finalNames = finalDimNames[0];

	int nMarkerPatternIDs = args.markerPatternData.allMarkerPatterns.size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations.begin(), args.intercrossingGenerations.end());
	int minSelfing = *std::min_element(args.selfingGenerations.begin(), args.selfingGenerations.end());
	int maxSelfing = *std::max_element(args.selfingGenerations.begin(), args.selfingGenerations.end());

	//This is basically just a huge lookup table
	allMarkerPairData<maxAlleles> computedContributions(nMarkerPatternIDs);
	
	constructLookupTableArgs<maxAlleles> lookupArgs(computedContributions, args.markerPatternData);
	lookupArgs.recombinationFractions = &args.recombinationFractions;
	lookupArgs.funnelEncodings = &args.funnelEncodings;
	lookupArgs.intercrossingGenerations = &args.intercrossingGenerations;
	lookupArgs.selfingGenerations = &args.selfingGenerations;
	constructLookupTable<nFounders, maxAlleles, infiniteSelfing>(lookupArgs);

#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static, 1)
#endif
	//This set of loops DOES grow with problem size. 
	for(int markerCounter1 = args.marker1Start; markerCounter1 < args.marker1End; markerCounter1++)
	{
		int markerPatternID1 = args.markerPatternData.markerPatternIDs[markerCounter1];
		for(int markerCounter2 = args.marker2Start; markerCounter2 < args.marker2End; markerCounter2++)
		{
			//For some bits in the lower triangle we have already calculated a corresponding bit in the upper triangle, so don't recalculate these. This means we have some values in args.result which are never set, but later code is aware of this
			if(markerCounter2 >= maxStart && markerCounter2 < minEnd && markerCounter1 >= maxStart && markerCounter1 < minEnd && markerCounter2 < markerCounter1) continue;
			int markerPatternID2 = args.markerPatternData.markerPatternIDs[markerCounter2];

			singleMarkerPairData<maxAlleles>& markerPairData = computedContributions(markerPatternID1, markerPatternID2);
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
				{
					int marker1Value = args.finals(finalCounter, markerCounter1);
					int marker2Value = args.finals(finalCounter, markerCounter2);
					if(marker1Value != NA_INTEGER && marker2Value != NA_INTEGER)
					{
						double contribution = 0;
						bool allowable = false;
						int intercrossingGenerations = args.intercrossingGenerations[finalCounter];
						int selfingGenerations = args.selfingGenerations[finalCounter];
						if(intercrossingGenerations == 0)
						{
							funnelID currentLineFunnelID = args.funnelIDs[finalCounter];
							allowable = markerPairData.allowableFunnel(currentLineFunnelID, selfingGenerations - minSelfing);
							if(allowable)
							{
								arrayType<maxAlleles>& perMarkerGenotypeValues = markerPairData.perFunnelData(recombCounter, currentLineFunnelID, 0);
								contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
							}
						}
						else if(intercrossingGenerations > 0)
						{
							allowable = markerPairData.allowableAI(intercrossingGenerations-1, selfingGenerations - minSelfing);
							if(allowable)
							{
								arrayType<maxAlleles>& perMarkerGenotypeValues = markerPairData.perAIGenerationData(recombCounter, intercrossingGenerations-1, 0);
								contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
							}
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
			std::fill(args.selfingGenerations.begin(), args.selfingGenerations.end(), 0);
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
/*		case 4:
			return estimateRfSpecificDesignInternal2<nFounders, 4>(args);
		case 5:
			return estimateRfSpecificDesignInternal2<nFounders, 5>(args);
		case 6:
			return estimateRfSpecificDesignInternal2<nFounders, 6>(args);
		case 7:
			return estimateRfSpecificDesignInternal2<nFounders, 7>(args);
		case 8:
			return estimateRfSpecificDesignInternal2<nFounders, 8>(args);*/
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
	bool hasAIC = *std::max_element(intercrossingGenerations.begin(), intercrossingGenerations.end()) > 0;

	//Check that everything has a unique funnel - For the case of the lines which are just selfing, we just check that one funnel. For AIC lines, we check the funnels of all the parent lines
	std::vector<funnelType> allFunnels;
	{
		std::vector<std::string> funnelWarnings, funnelErrors;
		estimateRfCheckFunnels(args.finals, args.founders, args.pedigree, intercrossingGenerations, funnelWarnings, funnelErrors, allFunnels);
		for(std::size_t errorIndex = 0; errorIndex < funnelErrors.size() && errorIndex < 6; errorIndex++)
		{
			Rprintf(funnelErrors[errorIndex].c_str());
		}
		if(funnelErrors.size() > 0)
		{
			return false;
		}
		for(std::size_t warningIndex = 0; warningIndex < funnelWarnings.size() && warningIndex < 6; warningIndex++)
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
	//We need to assign a unique ID to each marker pattern - Where by pattern we mean the combination of hetData and founder alleles. Note that this is possible because we just recoded everything to a standardised format.
	//Marker IDs are guaranteed to be contiguous numbers starting from 0 - So the set of all valid [0, markerPatterns.size()]. 
	//Note that markerPatternID is defined in unitTypes.hpp. It's just an integer (and automatically convertible to an integer), but it's represented by a unique type - This stops us from confusing it with an ordinary integer.
	markerPatternsToUniqueValuesArgs markerPatternConversionArgs;
	markerPatternConversionArgs.nFounders = nFounders;
	markerPatternConversionArgs.nMarkers = nMarkers;
	markerPatternConversionArgs.recodedFounders = recodedFounders;
	markerPatternConversionArgs.recodedHetData = recodedHetData;
	markerPatternsToUniqueValues(markerPatternConversionArgs);
	
	//map containing encodings of the funnels involved in the experiment (as key), and an associated unique index (again, using the encoded values directly is no good because they'll be all over the place). Unique indices are contiguous again.
	std::map<funnelEncoding, funnelID> funnelTranslation;
	//vector giving the funnel ID for each individual
	std::vector<funnelID> funnelIDs;
	//vector giving the encoded value for each individual
	std::vector<funnelEncoding> funnelEncodings;
	funnelsToUniqueValues(funnelTranslation, funnelIDs, funnelEncodings, allFunnels, nFounders);
	
	rfhaps_internal_args internal_args(args.lineWeights, args.recombinationFractions);
	internal_args.finals = recodedFinals;
	internal_args.pedigree = args.pedigree;
	internal_args.intercrossingGenerations.swap(intercrossingGenerations);
	internal_args.selfingGenerations.swap(selfingGenerations);
	internal_args.markerPatternData.swap(markerPatternConversionArgs);
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
