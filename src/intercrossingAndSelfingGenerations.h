#ifndef INTERCROSSING_AND_SELFING_GENERATIONS_HEADER
#define INTERCROSSING_AND_SELFING_GENERATIONS_HEADER
#include <Rcpp.h>
#include <vector>
//Work out the number of intercrossing generations for each individual in the final population
SEXP getIntercrossingAndSelfingGenerationsExport(SEXP pedigree_sexp, SEXP finals_sexp);
bool getIntercrossingAndSelfingGenerations(Rcpp::S4 pedigree, Rcpp::IntegerMatrix finals, int nFounders, std::vector<int>& intecrossing, std::vector<int>& selfing);
#endif
