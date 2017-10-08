#ifndef MATRIX_CHUNKS_HEADER_GUARD
#define MATRIX_CHUNKS_HEADER_GUARD
#include <Rcpp.h>
#include <vector>
class triangularIterator
{
public:
	triangularIterator(const std::vector<int>& markerRows, const std::vector<int>& markerColumns);
	std::pair<int, int> get() const;
	void next();
	bool isDone() const;
	triangularIterator& operator=(const triangularIterator& other);
private:
	const std::vector<int>& markerRows;
	const std::vector<int>& markerColumns;
	std::vector<int>::const_iterator markerRow, markerColumn;
};
class triangularIteratorPredicates
{
public:
	triangularIteratorPredicates(const std::vector<int>& markerRows, const std::vector<int>& markerColumns, std::function<bool(int)> predicate, std::function<bool(int, int)> jointPredicate);
	std::pair<int, int> get() const;
	void next();
	bool isDone() const;
	triangularIteratorPredicates& operator=(const triangularIteratorPredicates& other);
	unsigned long long getFlatIndex() const;
private:
	const std::vector<int>& markerRows;
	const std::vector<int>& markerColumns;
	std::vector<int>::const_iterator markerRow, markerColumn;
	std::function<bool(int)> predicate;
	std::function<bool(int, int)> jointPredicate;
};
SEXP countValuesToEstimateExported(SEXP markerRows, SEXP markerColumns);
unsigned long long countValuesToEstimate(const std::vector<int>& markerRows, const std::vector<int>& markerColumns);
SEXP singleIndexToPairExported(SEXP markerRows, SEXP markerColumns, SEXP index);
#endif
