#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "checkHets.h"
#include "generateGenotypes.h"
#include "alleleDataErrors.h"
extern "C"
{
	R_CallMethodDef callMethods[] = 
	{
		{"checkHets", (DL_FUNC)&checkHets, 1},
		{"generateGenotypes", (DL_FUNC)&generateGenotypes, 3},
		{"alleleDataErrors", (DL_FUNC)&alleleDataErrors, 2},
		{"listCodingErrors", (DL_FUNC)&listCodingErrors, 3},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap2(DllInfo *info)
	{
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}