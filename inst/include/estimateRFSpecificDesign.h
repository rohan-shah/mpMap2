#ifndef _ESTIMATE_RF_SPECIFIC_DESIGN_H
#define _ESTIMATE_RF_SPECIFIC_DESIGN_H
#include <Rcpp.h>
#include <vector>
struct estimateRFSpecificDesignArgs
{
	estimateRFSpecificDesignArgs(std::vector<double>& lineWeights, std::vector<double>& recombinationFractions)
		: recombinationFractions(recombinationFractions), lineWeights(lineWeights)
	{}
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 pedigree;
	Rcpp::S4 hetData;
	std::vector<double>& recombinationFractions;
	
	std::vector<double>& lineWeights;
	int marker1Start, marker1End;
	int marker2Start, marker2End;
	double* result;
	std::string error;
};
bool estimateRFSpecificDesign(estimateRFSpecificDesignArgs& args);
#endif

