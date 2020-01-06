#include "combineGenotypes.h"
SEXP combineGenotypes(SEXP Rfinals, SEXP RhetData)
{
	BEGIN_RCPP
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
		int nLines = finals.nrow();
		int nMarkers = finals.ncol()/2;
		Rcpp::IntegerMatrix result(nLines, nMarkers);
		for(Rcpp::List::iterator markerIterator = hetData.begin(); markerIterator != hetData.end(); markerIterator++)
		{
			Rcpp::IntegerMatrix hetDataThisMarker = *markerIterator;
			int markerIndex = std::distance(hetData.begin(), markerIterator);
			int hetDataRows = hetDataThisMarker.nrow();
			for(int lineCounter = 0; lineCounter < nLines; lineCounter++)
			{
				int allele1 = finals(lineCounter, markerIndex), allele2 = finals(lineCounter, markerIndex+nMarkers);
				int i = 0;
				for(; i < hetDataRows; i++)
				{
					if(hetDataThisMarker(i, 0) == allele1 && hetDataThisMarker(i, 1) == allele2)
					{
						result(lineCounter, markerIndex) = hetDataThisMarker(i, 2);
						break;
					}
				}
				if(i == hetDataRows)
				{
					throw std::runtime_error("Invalid marker encoding found. Please validate the input object using validObject.");
				}
			}
		}
		return Rcpp::wrap(result);
	END_RCPP
}

