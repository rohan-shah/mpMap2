#include "singleLocusProbabilities_RInterface.h"
#include "probabilities.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
template<int nFounders, bool infiniteSelfing> Rcpp::NumericVector singleLocusProbabilities(std::size_t nFunnels, int intercrossingGenerations, int selfingGenerations)
{
	array2<nFounders> data;
	if(intercrossingGenerations > 0)
	{
		singleLocusGenotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(data, selfingGenerations, nFunnels);
	}
	else
	{
		singleLocusGenotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(data, selfingGenerations, nFunnels);
	}
	Rcpp::NumericMatrix retVal(nFounders, nFounders);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			retVal(i, j) = data.values[i][j];
		}
	}
	return retVal;
}
SEXP singleLocusProbabilitiesFinite_RInterface(SEXP nFounders_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp, SEXP selfingGenerations_sexp)
{
BEGIN_RCPP
	int nFounders = Rcpp::as<int>(nFounders_sexp);
	int nFunnels = Rcpp::as<int>(nFunnels_sexp);
	int intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	int selfingGenerations = Rcpp::as<int>(selfingGenerations_sexp);
	if(nFounders == 2)
	{
		return singleLocusProbabilities<2, false>(nFunnels, intercrossingGenerations, selfingGenerations);
	}
	else if(nFounders == 4)
	{
		return singleLocusProbabilities<4, false>(nFunnels, intercrossingGenerations, selfingGenerations);
	}
	else if(nFounders == 8)
	{
		return singleLocusProbabilities<8, false>(nFunnels, intercrossingGenerations, selfingGenerations);
	}
	else if(nFounders == 16)
	{
		return singleLocusProbabilities<16, false>(nFunnels, intercrossingGenerations, selfingGenerations);
	}
	else 
	{
		throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
	}
END_RCPP
}
SEXP singleLocusProbabilitiesInfinite_RInterface(SEXP nFounders_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp)
{
BEGIN_RCPP
	int nFounders = Rcpp::as<int>(nFounders_sexp);
	int nFunnels = Rcpp::as<int>(nFunnels_sexp);
	int intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	if(nFounders == 2)
	{
		return singleLocusProbabilities<2, true>(nFunnels, intercrossingGenerations, 0);
	}
	else if(nFounders == 4)
	{
		return singleLocusProbabilities<4, true>(nFunnels, intercrossingGenerations, 0);
	}
	else if(nFounders == 8)
	{
		return singleLocusProbabilities<8, true>(nFunnels, intercrossingGenerations, 0);
	}
	else if(nFounders == 16)
	{
		return singleLocusProbabilities<16, true>(nFunnels, intercrossingGenerations, 0);
	}
	else 
	{
		throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
	}
END_RCPP
}
