#ifndef RECODE_FOUNDERS_AND_FINALS_AND_HET_DATA_HEADER_GUARD
#define RECODE_FOUNDERS_AND_FINALS_AND_HET_DATA_HEADER_GUARD
#include <Rcpp.h>
struct recodeDataStruct
{
	Rcpp::IntegerMatrix recodedFounders;
	Rcpp::IntegerMatrix recodedFinals;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix founders;
	Rcpp::IntegerMatrix finals;
	Rcpp::S4 hetData;
	unsigned int maxAlleles;
};
void recodeFoundersFinalsHets(recodeDataStruct& inputs);
#endif

