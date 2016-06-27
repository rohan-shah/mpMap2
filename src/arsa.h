#ifndef MPMAP2_ARSA_HEADER_GUARD
#define MPMAP2_ARSA_HEADER_GUARD
#include "Rcpp.h"
#include <vector>
#include <functional>
SEXP arsaExportedR(SEXP n_, SEXP dist_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_);
struct arsaArgs
{
public:
	R_xlen_t n;
	double* dist;
	int nReps;
	double temperatureMin, cool;
	double effortMultiplier;
	std::vector<int> bestPermutationAllReps;
	std::function<void(unsigned long,unsigned long)> progressFunction;
};
void arsaExported(R_xlen_t n, double* dist, int nReps, double temperatureMin, double cool, double effortMultiplier, std::vector<int>& bestPermutationAllReps, std::function<void(unsigned long,unsigned long)> progressFunction);
void arsa(arsaArgs& args);
#endif

