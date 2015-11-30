#include "matrixChunks.h"
SEXP countValuesToEstimateExported(SEXP marker1Start_, SEXP marker1End_, SEXP marker2Start_, SEXP marker2End_)
{
BEGIN_RCPP
	int marker1Start = Rcpp::as<int>(marker1Start_) - 1;
	int marker1End = Rcpp::as<int>(marker1End_) - 1;
	int marker2Start = Rcpp::as<int>(marker2Start_) - 1;
	int marker2End = Rcpp::as<int>(marker2End_) - 1;
	return Rcpp::wrap<int>(countValuesToEstimate(marker1Start, marker1End, marker2Start, marker2End));
END_RCPP
}
unsigned long long countValuesToEstimate(int marker1Start, int marker1End, int marker2Start, int marker2End)
{
	int marker1RangeSize = marker1End - marker1Start;
	int marker2RangeSize = marker2End - marker2Start;
	//This is true if the region of the upper triangle is in fact a rectangle (Doesn't intersect the diagonal)
	bool rectangularRegion = marker2Start >= marker1End-1;
	//Otherwise if this is true, the first row of the target region intersects the diagonal
	bool firstRowIntersects = marker2Start <= marker1Start;
	unsigned long long nValuesToEstimate;
	if(rectangularRegion) nValuesToEstimate = marker1RangeSize * marker2RangeSize;
	else
	{
		if(firstRowIntersects)
		{
			unsigned long long firstRowLength = marker2End - marker1Start;
			nValuesToEstimate = (firstRowLength * (firstRowLength + 1))/2;
			if(marker1End < marker2End) nValuesToEstimate -= ((marker2End - marker1End + 1) * (marker2End - marker1End))/2;
		}
		else
		{
			unsigned long long firstRowLength = marker2End - marker2Start;
			nValuesToEstimate = firstRowLength * (marker2Start - marker1Start) + (firstRowLength * (firstRowLength + 1))/2;
			if(marker1End < marker2End) nValuesToEstimate -= ((marker2End - marker1End) * (marker2End - marker1End+1))/2;
		}
	}
	return nValuesToEstimate;
}
SEXP singleIndexToPairExported(SEXP marker1Start_, SEXP marker1End_, SEXP marker2Start_, SEXP marker2End_, SEXP index_)
{
BEGIN_RCPP
	int marker1Start = Rcpp::as<int>(marker1Start_) - 1;
	int marker1End = Rcpp::as<int>(marker1End_) - 1;
	int marker2Start = Rcpp::as<int>(marker2Start_) - 1;
	int marker2End = Rcpp::as<int>(marker2End_) - 1;
	unsigned long long index = (unsigned long long)Rcpp::as<int>(index_) - 1;
	int markerCounter1, markerCounter2;
	singleIndexToPair(marker1Start, marker1End, marker2Start, marker2End, index, markerCounter1, markerCounter2);
	return Rcpp::IntegerVector::create(markerCounter2+1, markerCounter1+1);
END_RCPP
}
void singleIndexToPair(int marker1Start, int marker1End, int marker2Start, int marker2End, unsigned long long index, int& markerCounter1, int& markerCounter2)
{
	if(marker1End > marker2End) marker1End = marker2End;
	if(marker2Start < marker1Start) marker2Start = marker1Start;
	int marker1RangeSize = marker1End - marker1Start;
	int marker2RangeSize = marker2End - marker2Start;
	//This is true if the region of the upper triangle is in fact a rectangle (Doesn't intersect the diagonal)
	bool rectangularRegion = marker2Start >= marker1End-1;
	//Otherwise if this is true, the last column of the target region intersects the diagonal
	bool lastColumnIntersects = marker1End >= marker2End;
	if(rectangularRegion)
	{
		markerCounter1 = index % marker1RangeSize + marker1Start;
		markerCounter2 = index / marker1RangeSize + marker2Start;
	}
	else
	{
nonRectangular:
		unsigned long firstMissingColumnLength = 0;
		if(marker2Start >= marker1Start) firstMissingColumnLength = marker2Start - marker1Start;
		if(lastColumnIntersects)
		{
			unsigned long long lastColumnLength = marker2End - marker1Start;
			unsigned long long missingTriangular = (firstMissingColumnLength * (firstMissingColumnLength + 1))/2;

			//Zero-based column from the left side of the triangular region
			markerCounter2 = (std::size_t)sqrt(2 * (index + missingTriangular));
			if(((markerCounter2+1) * (markerCounter2 + 2))/2 < index + missingTriangular) throw std::runtime_error("Internal error");
			while((markerCounter2 * (markerCounter2 + 1)) / 2 > index + missingTriangular)
			{
				markerCounter2--;
				if(markerCounter2 < 0) throw std::runtime_error("Internal error");
			}
			markerCounter1 = index + missingTriangular - (markerCounter2 * (markerCounter2 + 1)) / 2;
			markerCounter1 += marker1Start;
			markerCounter2 += marker2Start - firstMissingColumnLength;
		}
		else
		{
			unsigned long long lastTriangularColumnLength = marker1End - marker1Start;
			unsigned long long missingTriangular = (firstMissingColumnLength * (firstMissingColumnLength + 1))/2;
			if(index + missingTriangular < (lastTriangularColumnLength * (lastTriangularColumnLength+1))/2)
			{
				lastColumnIntersects = true;
				goto nonRectangular;
			}
			index -= (lastTriangularColumnLength * (lastTriangularColumnLength+1))/2 - missingTriangular;
			markerCounter1 = index % marker1RangeSize + marker1Start;
			markerCounter2 = index / marker1RangeSize + marker2Start + lastTriangularColumnLength - firstMissingColumnLength;
		}
	}
}
