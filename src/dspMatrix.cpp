#include "dspMatrix.h"
#include "matrixChunks.h"
SEXP assignDspMatrixFromEstimateRF(SEXP destination_, SEXP rowIndices_, SEXP columnIndices_, SEXP source_)
{
BEGIN_RCPP
	Rcpp::S4 destination = destination_;
	Rcpp::NumericVector source = source_;
	Rcpp::NumericVector destinationData = destination.slot("x");
	Rcpp::IntegerVector rowIndices = rowIndices_;
	Rcpp::IntegerVector columnIndices = columnIndices_;

	if(&(source(0)) == &(destinationData(0)))
	{
		throw std::runtime_error("Source and destination cannot be the same in assignRawSymmetricMatrixDiagonal");
	}

	std::vector<int> markerRows, markerColumns;
	markerRows = Rcpp::as<std::vector<int> >(rowIndices);
	markerColumns = Rcpp::as<std::vector<int> >(columnIndices);
	if(countValuesToEstimate(markerRows, markerColumns) != (unsigned long long)source.size())
	{
		throw std::runtime_error("Mismatch between index length and source object size");
	}

	triangularIterator iterator(markerRows, markerColumns);
	int counter = 0;
	for(; !iterator.isDone(); iterator.next())
	{
		std::pair<int, int> markerPair = iterator.get();
		int markerRow = markerPair.first, markerColumn = markerPair.second;
		destinationData((markerColumn*(markerColumn-1))/2 + (markerRow - 1)) = source(counter);
		counter++;
	}
	return R_NilValue;
END_RCPP
}
