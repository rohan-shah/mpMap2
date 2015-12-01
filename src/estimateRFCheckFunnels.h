#ifndef ESTIMATE_RF_CHECK_FUNNELS_HEADER_GUARD
#define ESTIMATE_RF_CHECK_FUNNELS_HEADER_GUARD
#include <Rcpp.h>
#include <vector>
#include <string>
#include "estimateRF.h"
void estimateRFCheckFunnels(Rcpp::IntegerMatrix finals, Rcpp::IntegerMatrix founders, Rcpp::List hetData, Rcpp::S4 pedigree, std::vector<int>& intercrossingGenerations, std::vector<std::string>& warnings, std::vector<std::string>& errors, std::vector<funnelType>& allFunnels);
#endif