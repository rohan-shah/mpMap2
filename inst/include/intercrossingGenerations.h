#ifndef INTERCROSSING_GENERATIONS_HEADER
#define INTERCROSSING_GENERATIONS_HEADER
#include <Rcpp.h>
#include <vector>
//Work out the number of intercrossing generations for each individual in the final population
bool getIntercrossingGenerations(Rcpp::S4 pedigree, Rcpp::IntegerMatrix finals, int nFounders, std::vector<int>& output);
#endif
