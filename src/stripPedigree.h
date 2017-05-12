#ifndef STRIP_PEDIGREE_HEADER_GUARD_MPMAP2
#define STRIP_PEDIGREE_HEADER_GUARD_MPMAP2
#include <Rcpp.h>
#ifdef USE_BOOST
SEXP stripPedigree(SEXP pedigree_sexp, SEXP finalLines_sexp);
#endif
#endif
