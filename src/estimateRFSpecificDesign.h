#ifndef _ESTIMATE_RF_SPECIFIC_DESIGN_H
#define _ESTIMATE_RF_SPECIFIC_DESIGN_H
#include <Rcpp.h>
#include <vector>
#include "markerPatternsToUniqueValues.h"
#include "funnelsToUniqueValues.h"
#include "matrixChunks.h"
#include <functional>
struct estimateRFSpecificDesignArgs
{
	estimateRFSpecificDesignArgs(std::vector<double>& recombinationFractions)
		: recombinationFractions(recombinationFractions)
	{}
	estimateRFSpecificDesignArgs(const estimateRFSpecificDesignArgs& other)
		:founders(other.founders), finals(other.finals), pedigree(other.pedigree), hetData(other.hetData), recombinationFractions(other.recombinationFractions), lineWeights(other.lineWeights)
	{}
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	Rcpp::S4 hetData;
	std::vector<double>& recombinationFractions;
	
	std::vector<double> lineWeights;
};
struct rfhaps_internal_args
{
	rfhaps_internal_args(const std::vector<double>& recombinationFractions, triangularIterator& startPosition)
	: recombinationFractions(recombinationFractions), startPosition(startPosition)
	{}
	rfhaps_internal_args(rfhaps_internal_args&& other)
		:finals(other.finals), founders(other.founders), pedigree(other.pedigree), recombinationFractions(other.recombinationFractions), intercrossingGenerations(std::move(other.intercrossingGenerations)), selfingGenerations(std::move(other.selfingGenerations)), lineWeights(std::move(other.lineWeights)), markerPatternData(std::move(other.markerPatternData)), hasAI(other.hasAI), maxAlleles(other.maxAlleles), result(other.result), funnelIDs(std::move(other.funnelIDs)), funnelEncodings(std::move(other.funnelEncodings)), startPosition(other.startPosition)
	{}
	Rcpp::IntegerMatrix finals, founders;
	Rcpp::S4 pedigree;
	const std::vector<double>& recombinationFractions;
	
	std::vector<int> intercrossingGenerations;
	std::vector<int> selfingGenerations;
	std::vector<double> lineWeights;
	markerPatternsToUniqueValuesArgs markerPatternData;
	bool hasAI;
	//maximum number of marker alleles present
	int maxAlleles;
	//zeroed by the calling code. Must be added to, not overwritten. 
	double* result;
	unsigned long long valuesToEstimateInChunk;
	std::vector<funnelID> funnelIDs;
	std::vector<funnelEncoding> funnelEncodings;
	triangularIterator startPosition;
	std::function<void(unsigned long long)> updateProgress;
};
unsigned long long estimateLookup(rfhaps_internal_args& internal_args);
bool estimateRFSpecificDesign(rfhaps_internal_args& internal_args);
/* Preprocess inputs
 *
 * Preprocess inputs in preparation for estimating recombination fractions
 *
 * This function performs the following steps:
 * 1. Get out the number of generations of intercrossing and selfing, for each line. If we are assuming infinite generations of selfing, a value of 0 is used the selfing generations. This will be ignored later on. 
 * 2. Check for coding errors in the founders, finals and hetData. 
 * 3. Check that the data is consistent with the input funnels.
 * 4. Recode the founders and finals data in a normalised form, also working out the maximum number of alleles per marker, across the dataset
 * 5. Assign a unique value to each marker pattern and funnel. 
 * 6. In the case of infinite generations of selfing, replace hetrozygote values with NA. 
 * @param args Input arguments
 * @param internalArgs preprocessed arguments
 * @return A boolean value, with true indicating success. False indicates an error.
 */
bool toInternalArgs(estimateRFSpecificDesignArgs&& args, rfhaps_internal_args& internalArgs, std::string& error);
#endif

