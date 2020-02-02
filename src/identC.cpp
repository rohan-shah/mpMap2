#include "Rcpp.h"
//Copied from the methods package
Rcpp::RObject identC_(SEXP e1, SEXP e2)
{
    if(TYPEOF(e1) == STRSXP && TYPEOF(e2) == STRSXP &&
       LENGTH(e1) == 1 && LENGTH(e2) == 1 &&
       STRING_ELT(e1, 0) == STRING_ELT(e2, 0))
	return Rcpp::wrap(true);
    else
	return Rcpp::wrap(false);
}
SEXP identC(SEXP e1, SEXP e2)
{
	return identC_(e1, e2);
}

