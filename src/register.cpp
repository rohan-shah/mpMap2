#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "checkHets.h"
#include "generateGenotypes.h"
#include "alleleDataErrors.h"
#include "estimateRf.h"
extern "C"
{
	char* package_name = "mpMap2";
	void R_init_Rcpp(DllInfo *info);
	R_CallMethodDef callMethods[] = 
	{
		{"checkHets", (DL_FUNC)&checkHets, 1},
		{"generateGenotypes", (DL_FUNC)&generateGenotypes, 3},
		{"alleleDataErrors", (DL_FUNC)&alleleDataErrors, 2},
		{"listCodingErrors", (DL_FUNC)&listCodingErrors, 3},
		{"estimateRf", (DL_FUNC)&estimateRf, 7},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap2(DllInfo *info)
	{
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
		R_init_Rcpp(info);
	}
}