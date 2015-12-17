#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "checkHets.h"
#include "generateGenotypes.h"
#include "alleleDataErrors.h"
#include "estimateRF.h"
#include "internal.h"
#include "fourParentPedigreeRandomFunnels.h"
#include "matrixChunks.h"
#include "rawSymmetricMatrix.h"
#include "dspMatrix.h"
#include "preClusterStep.h"
#include "hclustMatrices.h"
extern "C"
{
	char* package_name = "mpMap2";
	R_CallMethodDef callMethods[] = 
	{
		{"checkHets", (DL_FUNC)&checkHets, 1},
		{"generateGenotypes", (DL_FUNC)&generateGenotypes, 3},
		{"alleleDataErrors", (DL_FUNC)&alleleDataErrors, 2},
		{"listCodingErrors", (DL_FUNC)&listCodingErrors, 3},
		{"estimateRF", (DL_FUNC)&estimateRF, 8},
		{"fourParentPedigreeRandomFunnels", (DL_FUNC)&fourParentPedigreeRandomFunnels, 4},
		{"countValuesToEstimate", (DL_FUNC)&countValuesToEstimateExported, 2},
		{"singleIndexToPair", (DL_FUNC)&singleIndexToPairExported, 3},
		{"rawSymmetricMatrixSubsetIndices", (DL_FUNC)&rawSymmetricMatrixSubsetIndices, 4},
		{"rawSymmetricMatrixSubsetObject", (DL_FUNC)&rawSymmetricMatrixSubsetObject, 2},
		{"assignRawSymmetricMatrixFromEstimateRF", (DL_FUNC)&assignRawSymmetricMatrixFromEstimateRF, 4},
		{"assignRawSymmetricMatrixDiagonal", (DL_FUNC)&assignRawSymmetricMatrixDiagonal, 3},
		{"assignDspMatrixFromEstimateRF", (DL_FUNC)&assignDspMatrixFromEstimateRF, 4},
		{"preClusterStep", (DL_FUNC)&preClusterStep, 1},
		{"hclustThetaMatrix", (DL_FUNC)&hclustThetaMatrix, 2},
		{"hclustCombinedMatrix", (DL_FUNC)&hclustCombinedMatrix, 2},
		{"hclustLodMatrix", (DL_FUNC)&hclustLodMatrix, 2},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap2(DllInfo *info)
	{
		std::vector<R_CallMethodDef> callMethodsVector;
		R_CallMethodDef* mpMap2CallMethods = callMethods;
		while(mpMap2CallMethods->name != NULL) mpMap2CallMethods++;
		callMethodsVector.insert(callMethodsVector.begin(), callMethods, mpMap2CallMethods);

#ifdef CUSTOM_STATIC_RCPP
		R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
		R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
		while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
		callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
#endif
		R_CallMethodDef blank = {NULL, NULL, 0};
		callMethodsVector.push_back(blank);

		R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
#ifdef CUSTOM_STATIC_RCPP
		init_Rcpp_cache();
#endif
	}
}
