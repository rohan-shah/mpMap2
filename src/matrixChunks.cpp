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
triangularIterator& triangularIterator::operator=(const triangularIterator& other)
{
	if(&markerRows != &other.markerRows || &markerColumns != &other.markerColumns) throw std::runtime_error("Internal error");
	markerRow = other.markerRow;
	markerColumn = other.markerColumn;
	return *this;
}
void triangularIterator::next()
{
	if(markerColumn == markerColumns.end()) throw std::runtime_error("Tried to increment iterator past the end");
	do
	{
		markerRow++;
		if(markerRow == markerRows.end())
		{
			markerRow = markerRows.begin();
			markerColumn++;
			if(markerColumn == markerColumns.end()) break;
		}
	}
	while(*markerRow > *markerColumn);
}
bool triangularIterator::isDone() const
{
	return markerColumn == markerColumns.end();
}
triangularIteratorPredicates::triangularIteratorPredicates(const std::vector<int>& markerRows, const std::vector<int>& markerColumns, std::function<bool(int)> predicate, std::function<bool(int, int)> jointPredicate)
	: markerRows(markerRows), markerColumns(markerColumns), markerRow(markerRows.begin()), markerColumn(markerColumns.begin()), predicate(predicate), jointPredicate(jointPredicate)
{
	if(*markerRow > *markerColumn || !jointPredicate(*markerRow, *markerColumn)) next();
}
std::pair<int, int> triangularIteratorPredicates::get() const
{
	return std::make_pair(*markerRow, *markerColumn);
}
triangularIteratorPredicates& triangularIteratorPredicates::operator=(const triangularIteratorPredicates& other)
{
	if(&markerRows != &other.markerRows || &markerColumns != &other.markerColumns) throw std::runtime_error("Internal error");
	markerRow = other.markerRow;
	markerColumn = other.markerColumn;
	predicate = other.predicate;
	jointPredicate = other.jointPredicate;
	return *this;
}
void triangularIteratorPredicates::next()
{
	if(markerColumn == markerColumns.end()) throw std::runtime_error("Tried to increment iterator past the end");
	do
	{
		markerRow++;
		if(markerRow == markerRows.end())
		{
			markerRow = markerRows.begin();
			markerColumn++;
			while(markerColumn != markerColumns.end() && !predicate(*markerColumn)) markerColumn++;
			if(markerColumn == markerColumns.end()) break;
		}
	}
	while(*markerRow > *markerColumn || !jointPredicate(*markerRow, *markerColumn));
}
bool triangularIteratorPredicates::isDone() const
{
	return markerColumn == markerColumns.end();
}
unsigned long long triangularIteratorPredicates::getFlatIndex() const
{
	unsigned long long column = std::distance(markerColumns.begin(), markerColumn);
	unsigned long long row = std::distance(markerRows.begin(), markerRow);
	return (column * (column + 1ULL) / 2ULL) + row;
}
SEXP countValuesToEstimateExported(SEXP markerRows_, SEXP markerColumns_)
{
BEGIN_RCPP
	std::vector<int> markerRows = Rcpp::as<std::vector<int> >(markerRows_);
	std::vector<int> markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
	unsigned long long result = countValuesToEstimate(markerRows, markerColumns);
	if (result > std::numeric_limits<int>::max()) throw std::runtime_error("Return value exceeds maximum value of 32-bit integer");
	return Rcpp::wrap<int>((int)result);
END_RCPP
}
unsigned long long countValuesToEstimate(const std::vector<int>& markerRows, const std::vector<int>& markerColumns)
{
	unsigned long long nValuesToEstimate = 0;
	std::vector<int> markerRowsCopied = markerRows, markerColumnsCopied = markerColumns;
	std::sort(markerRowsCopied.begin(), markerRowsCopied.end());
	std::sort(markerColumnsCopied.begin(), markerColumnsCopied.end());

	std::vector<int>::iterator columnIterator = markerColumnsCopied.begin();
	std::vector<int>::iterator rowIterator = markerRowsCopied.begin();
	while(columnIterator != markerColumnsCopied.end())
	{
		int column = *columnIterator;
		while(rowIterator != markerRowsCopied.end() && *rowIterator <= column)
		{
			rowIterator++;
		}
		nValuesToEstimate += std::distance(markerRowsCopied.begin(), rowIterator);
		columnIterator++;
	}
	return nValuesToEstimate;
}
SEXP singleIndexToPairExported(SEXP markerRows_, SEXP markerColumns_, SEXP index_)
{
BEGIN_RCPP
	std::vector<int> markerRows = Rcpp::as<std::vector<int> >(markerRows_);
	std::vector<int> markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
	triangularIterator iterator(markerRows, markerColumns);
	unsigned long long index = (unsigned long long)Rcpp::as<int>(index_) - 1;
	while(index > 0)
	{
		iterator.next();
		index--;
	}
	std::pair<int, int> markerPair = iterator.get();
	return Rcpp::IntegerVector::create(markerPair.first, markerPair.second);
END_RCPP
}

