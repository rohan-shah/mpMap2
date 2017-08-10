#ifndef ESTIMATE_RF_SINGLE_DESIGN_HEADER_GUARD
#define ESTIMATE_RF_SINGLE_DESIGN_HEADER_GUARD
#include <Rcpp.h>
SEXP estimateRFSingleDesign(SEXP object, SEXP recombinationFractions, SEXP markerRows, SEXP markerColumns, SEXP lineWeights, SEXP keepLod, SEXP keepLkhd, SEXP verbose);
#endif
