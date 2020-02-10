#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "checkHets.h"
#include "generateGenotypes.h"
#include "alleleDataErrors.h"
#include "estimateRF.h"
#ifdef CUSTOM_STATIC_RCPP
#include "internal.h"
#endif
#include "parsePurdy.h"
#include "fourParentPedigreeRandomFunnels.h"
#include "fourParentPedigreeSingleFunnel.h"
#include "eightParentPedigreeRandomFunnels.h"
#include "eightParentPedigreeSingleFunnel.h"
#include "eightParentPedigreeImproperFunnels.h"
#include "sixteenParentPedigreeRandomFunnels.h"
#include "matrixChunks.h"
#include "rawSymmetricMatrix.h"
#include "dspMatrix.h"
#include "preClusterStep.h"
#include "hclustMatrices.h"
#include "mpMap2_openmp.h"
#include "order.h"
#include "arsa.h"
#include "arsaRawREntryPoint.h"
#include "impute.h"
#include "multiparentSNP.h"
#include "imputeFounders.h"
#include "checkImputedBounds.h"
#include "generateDesignMatrix.h"
#include "compressedProbabilities_RInterface.h"
#include "testDistortion.h"
#include "removeHets.h"
#include "computeGenotypeProbabilities.h"
#include "computeAllEpistaticChiSquared.h"
#include "getAllFunnels.h"
#include "getAllFunnelsIncAIC.h"
#ifdef USE_BOOST
	#include "reorderPedigree.h"
#endif
#include "transformForMPWGAIM.h"
#include "intercrossingAndSelfingGenerations.h"
#include "stripPedigree.h"
#include "transposeProbabilities.h"
#include "assignFounderPattern.h"
#include "combineRFDisjoint.h"
#include "estimateRFSingleDesign.h"
#include "expandedProbabilities_RInterface.h"
#include "singleLocusProbabilities_RInterface.h"
#include "arsaRaw.h"
#include "identC.h"
extern "C"
{
	const char* package_name = "mpMap2";
	R_CallMethodDef callMethods[] = 
	{
		{"checkHets", (DL_FUNC)&checkHets, 1},
		{"generateGenotypes", (DL_FUNC)&generateGenotypes, 3},
		{"alleleDataErrors", (DL_FUNC)&alleleDataErrors, 2},
		{"listCodingErrors", (DL_FUNC)&listCodingErrors, 3},
		{"estimateRF", (DL_FUNC)&estimateRF, 9},
		{"fourParentPedigreeRandomFunnels", (DL_FUNC)&fourParentPedigreeRandomFunnels, 4},
		{"fourParentPedigreeSingleFunnel", (DL_FUNC)&fourParentPedigreeSingleFunnel, 4},
		{"eightParentPedigreeRandomFunnels", (DL_FUNC)&eightParentPedigreeRandomFunnels, 4},
		{"eightParentPedigreeSingleFunnel", (DL_FUNC)&eightParentPedigreeSingleFunnel, 4},
		{"sixteenParentPedigreeRandomFunnels", (DL_FUNC)&sixteenParentPedigreeRandomFunnels, 4},
		{"countValuesToEstimate", (DL_FUNC)&countValuesToEstimateExported, 2},
		{"singleIndexToPair", (DL_FUNC)&singleIndexToPairExported, 3},
		{"rawSymmetricMatrixSubsetIndices", (DL_FUNC)&rawSymmetricMatrixSubsetIndices, 4},
		{"rawSymmetricMatrixSubsetObject", (DL_FUNC)&rawSymmetricMatrixSubsetObject, 2},
		{"rawSymmetricMatrixToDist", (DL_FUNC)&rawSymmetricMatrixToDist, 1},
		{"constructDissimilarityMatrix", (DL_FUNC)&constructDissimilarityMatrix, 2},
		{"assignRawSymmetricMatrixFromEstimateRF", (DL_FUNC)&assignRawSymmetricMatrixFromEstimateRF, 4},
		{"assignRawSymmetricMatrixDiagonal", (DL_FUNC)&assignRawSymmetricMatrixDiagonal, 3},
		{"assignDspMatrixFromEstimateRF", (DL_FUNC)&assignDspMatrixFromEstimateRF, 4},
		{"preClusterStep", (DL_FUNC)&preClusterStep, 1},
		{"hclustThetaMatrix", (DL_FUNC)&hclustThetaMatrix, 2},
		{"hclustCombinedMatrix", (DL_FUNC)&hclustCombinedMatrix, 2},
		{"hclustLodMatrix", (DL_FUNC)&hclustLodMatrix, 2},
		{"omp_set_num_threads", (DL_FUNC)&mpMap2_omp_set_num_threads, 1},
		{"omp_get_num_threads", (DL_FUNC)&mpMap2_omp_get_num_threads, 0},
		{"order", (DL_FUNC)&order, 9},
		{"checkRawSymmetricMatrix", (DL_FUNC)&checkRawSymmetricMatrix, 1},
		{"arsa", (DL_FUNC)&arsaExportedR, 8},
		{"imputeGroup", (DL_FUNC)&imputeGroup, 4},
		{"multiparentSNPRemoveHets", (DL_FUNC)&multiparentSNPRemoveHets, 1},
		{"multiparentSNPKeepHets", (DL_FUNC)&multiparentSNPKeepHets, 1},
		{"rawSymmetricMatrixSubsetByMatrix", (DL_FUNC)&rawSymmetricMatrixSubsetByMatrix, 2},
		{"imputeFounders", (DL_FUNC)&imputeFounders, 7},
		{"checkImputedBounds", (DL_FUNC)&checkImputedBounds, 1},
		{"generateDesignMatrix", (DL_FUNC)&generateDesignMatrix, 2},
		{"compressedProbabilities", (DL_FUNC)&compressedProbabilities_RInterface, 6},
		{"eightParentPedigreeImproperFunnels", (DL_FUNC)&eightParentPedigreeImproperFunnels, 3},
#ifdef USE_BOOST
		{"reorderPedigree", (DL_FUNC)&reorderPedigree, 3},
		{"stripPedigree", (DL_FUNC)&stripPedigree, 2},
#endif
		{"testDistortion", (DL_FUNC)&testDistortion, 1},
		{"removeHets", (DL_FUNC)&removeHets, 3},
		{"computeGenotypeProbabilities", (DL_FUNC)&computeGenotypeProbabilities, 6},
		{"transformForMPWGAIM", (DL_FUNC)&transformForMPWGAIM, 1},
		{"parsePurdy", (DL_FUNC)&parsePurdy, 2}, 
		{"computeAllEpistaticChiSquared", (DL_FUNC)&computeAllEpistaticChiSquared, 4},
		{"getAllFunnels", (DL_FUNC)&getAllFunnels, 2}, 
		{"getAllFunnelsIncAIC", (DL_FUNC)&getAllFunnelsIncAIC, 2}, 
		{"rawSymmetricMatrixUncompress", (DL_FUNC)&rawSymmetricMatrixUncompress, 1},
		{"getIntercrossingAndSelfingGenerationsExport", (DL_FUNC)&getIntercrossingAndSelfingGenerationsExport, 2},
		{"transposeProbabilities", (DL_FUNC)&transposeProbabilities, 1},
		{"assignFounderPattern", (DL_FUNC)&assignFounderPattern, 2},
		{"combineRFDisjoint", (DL_FUNC)&combineRFDisjoint, 2},
		{"estimateRFSingleDesign", (DL_FUNC)&estimateRFSingleDesign, 8},
		{"expandedProbabilitiesInfinite", (DL_FUNC)&expandedProbabilitiesInfinite_RInterface, 4},
		{"expandedProbabilitiesFinite", (DL_FUNC)&expandedProbabilitiesFinite_RInterface, 6},
		{"singleLocusProbabilitiesFinite", (DL_FUNC)&singleLocusProbabilitiesFinite_RInterface, 4},
		{"singleLocusProbabilitiesInfinite", (DL_FUNC)&singleLocusProbabilitiesInfinite_RInterface, 4},
		{"identC", (DL_FUNC)&identC, 2},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap2(DllInfo *info)
	{
		R_useDynamicSymbols(info, FALSE);
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
		//Make non-const string
		char package_name_copy[100];
		memset(package_name_copy, '\0', sizeof(package_name_copy));
		strcpy(package_name_copy, package_name);
		R_RegisterCCallable(package_name_copy, "impute", (DL_FUNC)&impute);
		R_RegisterCCallable(package_name_copy, "constructDissimilarityMatrixInternal", (DL_FUNC)&constructDissimilarityMatrixInternal);
		R_RegisterCCallable(package_name_copy, "arsaRaw", (DL_FUNC)&arsaRaw::arsaRawExported);
		R_RegisterCCallable(package_name_copy, "arsa", (DL_FUNC)&arsa);
	}
}
