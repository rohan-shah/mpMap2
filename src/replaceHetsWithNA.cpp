#include "recodeHetsAsNA.h"
bool replaceHetsWithNA(Rcpp::IntegerMatrix recodedFounders, Rcpp::IntegerMatrix recodedFinals, Rcpp::List recodedHetData)
{
	int nFounders = recodedFounders.nrow();
	int nMarkers = recodedFounders.ncol();
	int nFinals = recodedFinals.nrow();
	bool retValue = false;
	std::vector<bool> isHet(100);
	for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		Rcpp::IntegerMatrix currentMarkerHetData = recodedHetData(markerCounter);
		int nAlleles = 0;
		for(int hetDataRowCounter = 0; hetDataRowCounter < currentMarkerHetData.nrow(); hetDataRowCounter++)
		{
			if(currentMarkerHetData(hetDataRowCounter, 0) != currentMarkerHetData(hetDataRowCounter, 1))
			{
				isHet[currentMarkerHetData(hetDataRowCounter, 2)] = true;
			}
			else
			{
				isHet[currentMarkerHetData(hetDataRowCounter, 2)] = false;
			}
		}
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if(isHet[recodedFinals(finalCounter, markerCounter)])
			{
				recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
				retValue = true;
			}
		}
	}
	return retValue;
}