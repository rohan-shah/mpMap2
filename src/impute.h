#ifndef IMPUTE_HEADER_GUARD_MPMAP2
#define IMPUTE_HEADER_GUARD_MPMAP2
#include <vector>
#include <string>
bool impute(unsigned char* theta, std::vector<double>& thetaLevels, double* lod, double* lkhd, std::vector<int>& markers, std::string& error);
#endif
