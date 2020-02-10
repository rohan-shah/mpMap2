#include "mpMap2_openmp.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
SEXP mpMap2_omp_set_num_threads(SEXP num)
{
BEGIN_RCPP
#ifdef USE_OPENMP
	int nThreads = Rcpp::as<int>(num);
	omp_set_num_threads(nThreads);
	return R_NilValue;
#else
	throw std::runtime_error("Not built with openmp suppport");
#endif
END_RCPP
}
SEXP mpMap2_omp_get_num_threads()
{
BEGIN_RCPP
#ifdef USE_OPENMP
	int nThreads = omp_get_num_threads();
	return Rcpp::wrap(nThreads);
#else
	throw std::runtime_error("Not built with openmp suppport");
#endif
END_RCPP
}
