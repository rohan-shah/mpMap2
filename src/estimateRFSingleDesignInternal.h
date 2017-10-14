#ifndef ESTIMATE_RF_SINGLE_DESIGN_INTERNAL
#define ESTIMATE_RF_SINGLE_DESIGN_INTERNAL
#include <vector>
#include "markerPatternsToUniqueValues.h"
#include "funnelsToUniqueValues.h"
#include <set>
struct estimateRFSingleDesignArgs
{
	estimateRFSingleDesignArgs(std::vector<double>& recombinationFractions)
		: recombinationFractions(recombinationFractions)
	{}
	estimateRFSingleDesignArgs(const estimateRFSingleDesignArgs& other)
		:founders(other.founders), finals(other.finals), pedigree(other.pedigree), hetData(other.hetData), recombinationFractions(other.recombinationFractions), lineWeights(other.lineWeights)
	{}
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	Rcpp::S4 hetData;
	std::vector<double>& recombinationFractions;
	std::vector<double> lineWeights;
};
struct estimateRFSingleDesignInternalArgs
{
	estimateRFSingleDesignInternalArgs(const std::vector<double>& recombinationFractions)
	: recombinationFractions(recombinationFractions), markerRows(NULL), markerColumns(NULL), lod(NULL), lkhd(NULL), theta(NULL)
	{}
	estimateRFSingleDesignInternalArgs(estimateRFSingleDesignInternalArgs&& other)
		:finals(other.finals), founders(other.founders), pedigree(other.pedigree), recombinationFractions(other.recombinationFractions), intercrossingGenerations(std::move(other.intercrossingGenerations)), selfingGenerations(std::move(other.selfingGenerations)), lineWeights(std::move(other.lineWeights)), markerPatternData(std::move(other.markerPatternData)), hasAI(other.hasAI), lineFunnelIDs(std::move(other.lineFunnelIDs)), lineFunnelEncodings(std::move(other.lineFunnelEncodings)), allFunnelEncodings(std::move(other.allFunnelEncodings)), rowPatterns(std::move(other.rowPatterns)), columnPatterns(std::move(other.columnPatterns)), markerRows(other.markerRows), markerColumns(other.markerColumns), lod(other.lod), lkhd(other.lkhd)
	{}
	Rcpp::IntegerMatrix finals, founders;
	Rcpp::S4 pedigree;
	const std::vector<double>& recombinationFractions;
	
	std::vector<int> intercrossingGenerations;
	std::vector<int> selfingGenerations;
	std::vector<double> lineWeights;
	markerPatternsToUniqueValuesArgs markerPatternData;
	bool hasAI;
	std::vector<funnelID> lineFunnelIDs;
	std::vector<funnelEncoding> lineFunnelEncodings;
	std::vector<funnelEncoding> allFunnelEncodings;
	std::function<void(unsigned long long)> updateProgress;
	std::vector<int> rowPatterns, columnPatterns;
	const std::vector<int>* markerRows, *markerColumns;
	double* lod, *lkhd;
	unsigned char* theta;
	std::vector<std::pair<int, int> > uniquePatternPairs;
};
bool toInternalArgs(estimateRFSingleDesignArgs&& args, estimateRFSingleDesignInternalArgs& internalArgs, std::string& error);
bool estimateRFSingleDesignInternal(estimateRFSingleDesignInternalArgs& args);
#endif
