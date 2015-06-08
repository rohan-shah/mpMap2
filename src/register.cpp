#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "checkHets.h"
#include "generateGenotypes.h"
#include "alleleDataErrors.h"
#include "estimateRF.h"
extern "C"
{
	char* package_name = "mpMap2";
#ifdef CUSTOM_STATIC_RCPP
	void R_init_Rcpp(DllInfo *info);
#endif
	R_CallMethodDef callMethods[] = 
	{
		{"checkHets", (DL_FUNC)&checkHets, 1},
		{"generateGenotypes", (DL_FUNC)&generateGenotypes, 3},
		{"alleleDataErrors", (DL_FUNC)&alleleDataErrors, 2},
		{"listCodingErrors", (DL_FUNC)&listCodingErrors, 3},
		{"estimateRF", (DL_FUNC)&estimateRF, 7},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap2(DllInfo *info)
	{
#ifdef CUSTOM_STATIC_RCPP
		R_init_Rcpp(info);
#endif
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}