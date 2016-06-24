#ifndef MPMAP2_ARSA_HEADER_GUARD
#define MPMAP2_ARSA_HEADER_GUARD
#include "Rcpp.h"
#include <vector>
SEXP arsaExportedR(SEXP n_, SEXP dist_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_);
struct arsaArgs
{
public:
	R_xlen_t n;
	double* dist;
	int nReps;
	double temperatureMin, cool;
	std::vector<int> bestPermutationAllReps;
};
void arsaExported(R_xlen_t n, double* dist, int nReps, double temperatureMin, double cool, std::vector<int>& bestPermutationAllReps);
void arsa(arsaArgs& args);
#endif

