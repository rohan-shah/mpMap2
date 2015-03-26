#include "combineGenotypes.h"
SEXP combineGenotypes(SEXP Rfinals, SEXP RhetData);
{
	RCPP_BEGIN
		Rcpp::IntegerMatrix finals;
		try
		{
			finals = Rfinals;
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input finals must be an integer matrix");
		}
		Rcpp::List hetData;
		try
		{
			hetData = RhetData;
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input hetData must be an S4 object");
		}
		int nLines = finals.nrows();
		int nMarkers = finals.ncols()/2;
		Rcpp::IntegerMatrix result(nLines, nMarkers);
		for(Rcpp::List::iterator markerIterator = hetData.begin(); markerIterator != hetData.end(); markerIterator++)
		{
			Rcpp::IntegerMatrix hetDataThisMarker = *markerIterator;
			int hetDataRows = hetDataThisMarker.nrows();
			for(int lineCounter = 0; lineCounter < nLines; lineCounter++)
			{
				int relevantRow = 0;
				int allele1 = 
			}
		}
	RCPP_END
}