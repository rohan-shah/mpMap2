#ifndef MATRIX_CHUNKS_HEADER_GUARD
#define MATRIX_CHUNKS_HEADER_GUARD
#include <Rcpp.h>
SEXP countValuesToEstimateExported(SEXP marker1Start, SEXP marker1End, SEXP marker2Start, SEXP marker2End);
unsigned long long countValuesToEstimate(int marker1Start, int marker1End, int marker2Start, int marker2End);
SEXP singleIndexToPairExported(SEXP marker1Start, SEXP marker1End, SEXP marker2Start, SEXP marker2End, SEXP index);
void singleIndexToPair(int marker1Start, int marker1End, int marker2Start, int marker2End, unsigned long long index, int& markerCounter1, int& markerCounter2);
#endif
