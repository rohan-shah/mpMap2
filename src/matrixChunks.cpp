#include "matrixChunks.h"
triangularIterator::triangularIterator(const std::vector<int>& markerRows, const std::vector<int>& markerColumns)
	: markerRows(markerRows), markerColumns(markerColumns), markerRow(markerRows.begin()), markerColumn(markerColumns.begin())
{
	if(*markerRow > *markerColumn) next();
}
std::pair<int, int> triangularIterator::get() const
{
	return std::make_pair(*markerRow, *markerColumn);
}
void triangularIterator::next()
{
	if(markerColumn == markerColumns.end()) throw std::runtime_error("Tried it increment iterator past the end");
	markerRow++;
	if(markerRow == markerRows.end() || *markerRow > *markerColumn)
	{
		markerRow = markerRows.begin();
		markerColumn++;
	}
}
bool triangularIterator::isDone() const
{
	return markerColumn == markerColumns.end();
}
SEXP countValuesToEstimateExported(SEXP markerRows_, SEXP markerColumns_)
{
BEGIN_RCPP
	std::vector<int> markerRows = Rcpp::as<std::vector<int> >(markerRows_);
	std::vector<int> markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
	for(std::vector<int>::iterator markerRow = markerRows.begin(); markerRow != markerRows.end(); markerRow++)
	{
		(*markerRow)--;
	}
	for(std::vector<int>::iterator markerColumn = markerColumns.begin(); markerColumn != markerColumns.end(); markerColumn++)
	{
		(*markerColumn)--;
	}
	std::sort(markerRows.begin(), markerRows.end());
	std::sort(markerColumns.begin(), markerColumns.end());
	return Rcpp::wrap<int>(countValuesToEstimate(markerRows, markerColumns));
END_RCPP
}
unsigned long long countValuesToEstimate(const std::vector<int>& markerRows, const std::vector<int>& markerColumns)
{
	unsigned long long nValuesToEstimate = 0;
	if(!std::is_sorted(markerRows.begin(), markerRows.end()) || !std::is_sorted(markerColumns.begin(), markerColumns.end()))
	{
		throw std::runtime_error("Inputs to countValuesToEstimate must be sorted");
	}
	std::vector<int>::const_iterator columnIterator = markerColumns.begin();
	std::vector<int>::const_iterator rowIterator = markerRows.begin();
	while(columnIterator != markerColumns.end())
	{
		int column = *columnIterator;
		while(rowIterator != markerRows.end() && *rowIterator <= column)
		{
			rowIterator++;
		}
		nValuesToEstimate += std::distance(markerRows.begin(), rowIterator);
		columnIterator++;
	}
	return nValuesToEstimate;
}
SEXP singleIndexToPairExported(SEXP markerRows_, SEXP markerColumns_, SEXP index_)
{
BEGIN_RCPP
	std::vector<int> markerRows = Rcpp::as<std::vector<int> >(markerRows_);
	std::vector<int> markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
	for(std::vector<int>::iterator markerRow = markerRows.begin(); markerRow != markerRows.end(); markerRow++)
	{
		(*markerRow)--;
	}
	for(std::vector<int>::iterator markerColumn = markerColumns.begin(); markerColumn != markerColumns.end(); markerColumn++)
	{
		(*markerColumn)--;
	}
	triangularIterator iterator(markerRows, markerColumns);
	unsigned long long index = (unsigned long long)Rcpp::as<int>(index_) - 1;
	while(index > 0)
	{
		iterator.next();
		index--;
	}
	std::pair<int, int> markerPair = iterator.get();
	return Rcpp::IntegerVector::create(markerPair.first+1, markerPair.second+1);
END_RCPP
}

