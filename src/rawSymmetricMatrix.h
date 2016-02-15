#ifndef RAW_SYMMETRIC_MATRIX_HEADER_GUARD
#define RAW_SYMMETRIC_MATRIX_HEADER_GUARD
#include <Rcpp.h>
SEXP rawSymmetricMatrixSubsetIndices(SEXP object, SEXP i, SEXP j, SEXP drop);
SEXP rawSymmetricMatrixSubsetObject(SEXP object, SEXP indices);
SEXP assignRawSymmetricMatrixFromEstimateRF(SEXP destination, SEXP rowIndices, SEXP columnIndices, SEXP source);
SEXP assignRawSymmetricMatrixDiagonal(SEXP destination, SEXP indices, SEXP source);
SEXP checkRawSymmetricMatrix(SEXP rawSymmetric);
SEXP rawSymmetricMatrixSubsetByMatrix(SEXP object_, SEXP index_);
#endif
