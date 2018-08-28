#ifndef IMPUTE_HEADER_GUARD_MPMAP2
#define IMPUTE_HEADER_GUARD_MPMAP2
#include <vector>
#include <string>
#include <functional>
#include <Rcpp.h>
bool impute(const unsigned char* originalTheta, unsigned char* imputedTheta, std::vector<double>& thetaLevels, double* lod, double* lkhd, std::vector<int>& markers, std::function<void(unsigned long, unsigned long)> statusFunction, bool allErrors, std::vector<std::pair<int, int> >& reportedErrors);
SEXP imputeGroup(SEXP mpcrossLG_sexp, SEXP verbose_sexp, SEXP group_sexp, SEXP allErrors);
#endif
