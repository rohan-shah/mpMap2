#ifndef RAW_SYMMETRIC_MATRIX_HEADER_GUARD
#define RAW_SYMMETRIC_MATRIX_HEADER_GUARD
#include <Rcpp.h>
SEXP rawSymmetricMatrixUncompress(SEXP object_);
SEXP rawSymmetricMatrixSubsetIndices(SEXP object, SEXP i, SEXP j, SEXP drop);
SEXP rawSymmetricMatrixSubsetObject(SEXP object, SEXP indices);
SEXP assignRawSymmetricMatrixFromEstimateRF(SEXP destination, SEXP rowIndices, SEXP columnIndices, SEXP source);
SEXP assignRawSymmetricMatrixDiagonal(SEXP destination, SEXP indices, SEXP source);
SEXP checkRawSymmetricMatrix(SEXP rawSymmetric);
SEXP rawSymmetricMatrixSubsetByMatrix(SEXP object_, SEXP index_);
SEXP rawSymmetricMatrixToDist(SEXP object);
SEXP constructDissimilarityMatrixInternal(unsigned char* data, std::vector<double>& levels, int size, SEXP clusters_, int start, const std::vector<int>& permutation);
SEXP constructDissimilarityMatrix(SEXP object, SEXP clusters);
#endif
