#ifndef REORDER_PEDIGREE_HEADER_GUARD
#define REORDER_PEDIGREE_HEADER_GUARD
#include <Rcpp.h>
#ifdef USE_BOOST
SEXP reorderPedigree(SEXP mother_sexp, SEXP father_sexp, SEXP lineNames_sexp);
#endif
#endif
