#ifndef _RFHAPS_H
#define _RFHAPS_H
#include <Rcpp.h>
struct funnelType
{
	int val[16];
};
SEXP estimateRF(SEXP object_, SEXP recombinationFractions_, SEXP marker1Range_, SEXP marker2Range_, SEXP lineWeights_, SEXP keepLod_, SEXP keepLkhd_);
#endif