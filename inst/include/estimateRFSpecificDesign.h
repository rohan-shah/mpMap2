#ifndef _ESTIMATE_RF_SPECIFIC_DESIGN_H
#define _ESTIMATE_RF_SPECIFIC_DESIGN_H
#include <Rcpp.h>
#include <vector>
#include "markerPatternsToUniqueValues.h"
#include "funnelsToUniqueValues.h"
struct estimateRFSpecificDesignArgs
{
	estimateRFSpecificDesignArgs(std::vector<double>& recombinationFractions, const std::vector<int>& markerRows, const std::vector<int>& markerColumns)
		: recombinationFractions(recombinationFractions), markerRows(markerRows), markerColumns(markerColumns)
	{}
	estimateRFSpecificDesignArgs(const estimateRFSpecificDesignArgs& other)
		:founders(other.founders), finals(other.finals), pedigree(other.pedigree), hetData(other.hetData), recombinationFractions(other.recombinationFractions), lineWeights(other.lineWeights), markerRows(other.markerRows), markerColumns(other.markerColumns), result(other.result), error(other.error)
	{}
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	Rcpp::S4 hetData;
	std::vector<double>& recombinationFractions;
	
	std::vector<double> lineWeights;
	const std::vector<int>& markerRows;
	const std::vector<int>& markerColumns;
	double* result;
	std::string error;
};
struct rfhaps_internal_args
{
	rfhaps_internal_args(std::vector<double>& lineWeights, std::vector<double>& recombinationFractions, const std::vector<int>& markerRows, const std::vector<int>& markerColumns)
	: recombinationFractions(recombinationFractions), lineWeights(lineWeights), markerRows(markerRows), markerColumns(markerColumns)
	{}
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	std::vector<double>& recombinationFractions;
	
	std::vector<int> intercrossingGenerations;
	std::vector<int> selfingGenerations;
	std::vector<double>& lineWeights;
	markerPatternsToUniqueValuesArgs markerPatternData;
	bool hasAI;
	const std::vector<int>& markerRows;
	const std::vector<int>& markerColumns;
	//maximum number of marker alleles present
	int maxAlleles;
	//zeroed by the calling code. Must be added to, not overwritten. 
	double* result;
	std::vector<funnelID> funnelIDs;
	std::vector<funnelEncoding> funnelEncodings;
};
unsigned long long estimateLookup(estimateRFSpecificDesignArgs& args);
bool estimateRFSpecificDesign(estimateRFSpecificDesignArgs& args);
bool toInternalArgs(estimateRFSpecificDesignArgs& args, rfhaps_internal_args& internalArgs, bool supressOutput);
#endif

