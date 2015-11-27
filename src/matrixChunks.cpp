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
			if(marker1End < marker2End) nValuesToEstimate -= ((marker1End - marker2End) * (marker1End - marker2End+1))/2;
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
	int marker1RangeSize = marker1End - marker1Start;
	int marker2RangeSize = marker2End - marker2Start;
	//This is true if the region of the upper triangle is in fact a rectangle (Doesn't intersect the diagonal)
	bool rectangularRegion = marker2Start >= marker1End-1;
	//Otherwise if this is true, the first row of the target region intersects the diagonal
	bool firstRowIntersects = marker2Start <= marker1Start;
	if(rectangularRegion)
	{
		markerCounter1 = index / marker2RangeSize + marker1Start;
		markerCounter2 = index % marker2RangeSize + marker2Start;
	}
	else
	{
nonRectangular:
		if(firstRowIntersects)
		{
			unsigned long long firstRowLength = marker2End - marker1Start;
			unsigned long long completeTriangular = (firstRowLength * (firstRowLength + 1))/2;
			unsigned long firstMissingRowLength = 0;
			if(marker1End < marker2End) firstMissingRowLength = marker2End - marker1End;

			//Zero-based index from the base of the triangular region
			int inversed = completeTriangular - index - 1;
			//Zero-based row above the base of the triangular region
			markerCounter1 = (std::size_t)sqrt(2 * inversed);
			if((markerCounter1 * (markerCounter1 + 1))/2 < inversed) throw std::runtime_error("Internal error");
			while((markerCounter1 * (markerCounter1 + 1)) / 2 > inversed)
			{
				markerCounter1--;
				if(markerCounter1 < 0) throw std::runtime_error("Internal error");
			}
			//Zero-based column index, but from the right-hand side (and we want from the left)
			markerCounter2 = inversed - (markerCounter1 * (markerCounter1 + 1)) / 2;

			markerCounter2 = firstMissingRowLength + markerCounter1 - markerCounter2;
			markerCounter1 = marker1RangeSize - markerCounter1 + marker1Start - 1;

			markerCounter2 += (markerCounter1 - marker1Start);
		}
		else
		{
			unsigned long long firstRowLength = marker2End - marker2Start;
			if(index > firstRowLength * (marker1Start - marker1Start))
			{
				//make some changes and jump back to the first case
				index -= firstRowLength * (marker1Start - marker1Start);
				firstRowIntersects = true;
				marker1Start = marker2Start;
				marker1End = marker2End;
				goto nonRectangular;
			}
			else
			{
				markerCounter1 = index / marker1RangeSize;
				markerCounter2 = index % marker1RangeSize;
			}
		}
	}
}
