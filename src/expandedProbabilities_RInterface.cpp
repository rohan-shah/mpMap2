#include "expandedProbabilities_RInterface.h"
#include "probabilities.hpp"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
template<int nFounders> Rcpp::NumericVector getCompressedProbabilities1(double r, int nFunnels, int intercrossingGenerations)
{
	std::array<double, compressedProbabilities<nFounders, true>::nDifferentProbs> result;
	const int selfingGenerations = 0;
	if(intercrossingGenerations > 0)
	{
		genotypeProbabilitiesWithIntercross<nFounders, true>(result, intercrossingGenerations, r, selfingGenerations, nFunnels);
	}
	else
	{
		genotypeProbabilitiesNoIntercross<nFounders, true>(result, r, selfingGenerations, nFunnels);
	}
	Rcpp::NumericMatrix retVal(nFounders, nFounders);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			retVal(i, j) = result[probabilityData<nFounders>::infiniteMask[i][j]];
		}
	}
	return retVal;
}
SEXP expandedProbabilitiesInfinite_RInterface(SEXP nFounders_sexp, SEXP r_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp)
{
BEGIN_RCPP
	int nFounders = Rcpp::as<int>(nFounders_sexp);
	double r = Rcpp::as<double>(r_sexp);
	int nFunnels = Rcpp::as<int>(nFunnels_sexp);
	int intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	if(nFounders == 2)
	{
		return getCompressedProbabilities1<2>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 4)
	{
		return getCompressedProbabilities1<4>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 8)
	{
		return getCompressedProbabilities1<8>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 16)
	{
		return getCompressedProbabilities1<16>(r, nFunnels, intercrossingGenerations);
	}
	else 
	{
		throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
	}
END_RCPP
}
