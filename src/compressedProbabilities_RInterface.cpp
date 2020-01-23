#include "compressedProbabilities_RInterface.h"
#include "probabilities.h"
template<int nFounders, bool infiniteSelfing> Rcpp::NumericVector getCompressedProbabilities2(double r, int nFunnels, int intercrossingGenerations, int selfingGenerations)
{
	std::array<double, compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs> result;
	if(intercrossingGenerations > 0)
	{
		genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(result, intercrossingGenerations, r, selfingGenerations, nFunnels);
	}
	else
	{
		genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(result, r, selfingGenerations, nFunnels);
	}
	Rcpp::NumericVector retVal(compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs);
	memcpy(&(retVal(0)), &(result[0]), sizeof(double)*compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs);
	return retVal;
}
template<int nFounders> Rcpp::NumericVector getCompressedProbabilities1(double r, int nFunnels, int intercrossingGenerations, int selfingGenerations, bool infiniteSelfing)
{
	if(infiniteSelfing)
	{
		return getCompressedProbabilities2<nFounders, true>(r, nFunnels, intercrossingGenerations, selfingGenerations);
	}
	return getCompressedProbabilities2<nFounders, false>(r, nFunnels, intercrossingGenerations, selfingGenerations);
}
SEXP compressedProbabilities_RInterface(SEXP nFounders_sexp, SEXP r_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp, SEXP selfingGenerations_sexp, SEXP infiniteSelfing_sexp)
{
BEGIN_RCPP
	int nFounders = Rcpp::as<int>(nFounders_sexp);
	double r = Rcpp::as<double>(r_sexp);
	int nFunnels = Rcpp::as<int>(nFunnels_sexp);
	int intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	int selfingGenerations = Rcpp::as<int>(selfingGenerations_sexp);
	bool infiniteSelfing = Rcpp::as<bool>(infiniteSelfing_sexp);
	if(nFounders == 2)
	{
		return getCompressedProbabilities1<2>(r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing);
	}
	else if(nFounders == 4)
	{
		return getCompressedProbabilities1<4>(r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing);
	}
	else if(nFounders == 8)
	{
		return getCompressedProbabilities1<8>(r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing);
	}
	else if(nFounders == 16)
	{
		return getCompressedProbabilities1<16>(r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing);
	}
	else 
	{
		throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
	}
END_RCPP
}
