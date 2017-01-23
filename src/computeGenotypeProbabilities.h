#ifndef COMPUTE_GENOTYPE_PROBABILITIES_HEADER_GUARD
#define COMPUTE_GENOTYPE_PROBABILITIES_HEADER_GUARD
#include "Rcpp.h"
SEXP computeGenotypeProbabilities(SEXP geneticData_sexp, SEXP map_sexp, SEXP homozygoteMissingProb_sexp, SEXP hetrozygoteMissingProb_sexp, SEXP errorProb_sexp, SEXP extraPositions_sexp);
#endif
