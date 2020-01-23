#include "expandedProbabilities_RInterface.h"
#include "probabilities.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
template<int nFounders> Rcpp::NumericVector getCompressedProbabilitiesInfinite(double r, int nFunnels, int intercrossingGenerations)
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
template<int nFounders> Rcpp::NumericVector getCompressedProbabilitiesFiniteUnphased(double r, int nFunnels, int intercrossingGenerations, int selfingGenerations)
{
	std::array<double, compressedProbabilities<nFounders, false>::nDifferentProbs> result;
	expandedProbabilitiesFiniteSelfing<nFounders> expandedProbabilities;
	if(intercrossingGenerations > 0)
	{
		genotypeProbabilitiesWithIntercross<nFounders, false>(result, intercrossingGenerations, r, selfingGenerations, nFunnels);
		expandedGenotypeProbabilitiesUnphased<nFounders, false, false>::withIntercross(expandedProbabilities, intercrossingGenerations, r, selfingGenerations, nFunnels);
	}
	else
	{
		genotypeProbabilitiesNoIntercross<nFounders, false>(result, r, selfingGenerations, nFunnels);
		expandedGenotypeProbabilitiesUnphased<nFounders, false, false>::noIntercross(expandedProbabilities, r, selfingGenerations, nFunnels);
	}
	Rcpp::NumericMatrix retVal(nFounders * nFounders, nFounders * nFounders);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			for(int k = 0; k < nFounders; k++)
			{
				for(int l = 0; l < nFounders; l++)
				{
					retVal[i * nFounders * nFounders * nFounders + j * nFounders * nFounders + k * nFounders + l] = expandedProbabilities.values[i][j][k][l];
				}
			}
		}
		
	}
	return retVal;
}
template<int nFounders> Rcpp::NumericVector getCompressedProbabilitiesFinitePhased(double r, int nFunnels, int intercrossingGenerations, int selfingGenerations)
{
	std::array<double, compressedProbabilities<nFounders, false>::nDifferentProbs> result;
	expandedProbabilitiesFiniteSelfing<nFounders> expandedProbabilities;
	if(intercrossingGenerations > 0)
	{
		genotypeProbabilitiesWithIntercross<nFounders, false>(result, intercrossingGenerations, r, selfingGenerations, nFunnels);
		expandedGenotypeProbabilitiesPhased<nFounders, false, false>::withIntercross(expandedProbabilities, intercrossingGenerations, r, selfingGenerations, nFunnels);
	}
	else
	{
		genotypeProbabilitiesNoIntercross<nFounders, false>(result, r, selfingGenerations, nFunnels);
		expandedGenotypeProbabilitiesPhased<nFounders, false, false>::noIntercross(expandedProbabilities, r, selfingGenerations, nFunnels);
	}
	Rcpp::NumericMatrix retVal(nFounders * nFounders, nFounders * nFounders);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			for(int k = 0; k < nFounders; k++)
			{
				for(int l = 0; l < nFounders; l++)
				{
					retVal[i * nFounders * nFounders * nFounders + j * nFounders * nFounders + k * nFounders + l] = expandedProbabilities.values[i][j][k][l];
				}
			}
		}
		
	}
	return retVal;
}
SEXP expandedProbabilitiesFinite_RInterface(SEXP nFounders_sexp, SEXP r_sexp, SEXP nFunnels_sexp, SEXP intercrossingGenerations_sexp, SEXP selfingGenerations_sexp, SEXP phased_sexp)
{
BEGIN_RCPP
	int nFounders = Rcpp::as<int>(nFounders_sexp);
	double r = Rcpp::as<double>(r_sexp);
	int nFunnels = Rcpp::as<int>(nFunnels_sexp);
	int intercrossingGenerations = Rcpp::as<int>(intercrossingGenerations_sexp);
	int selfingGenerations = Rcpp::as<int>(selfingGenerations_sexp);
	bool phased = Rcpp::as<bool>(phased_sexp);
	if(phased)
	{
		if(nFounders == 2)
		{
			return getCompressedProbabilitiesFinitePhased<2>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 4)
		{
			return getCompressedProbabilitiesFinitePhased<4>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 8)
		{
			return getCompressedProbabilitiesFinitePhased<8>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 16)
		{
			return getCompressedProbabilitiesFinitePhased<16>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else 
		{
			throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
		}
	}
	else
	{
		if(nFounders == 2)
		{
			return getCompressedProbabilitiesFiniteUnphased<2>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 4)
		{
			return getCompressedProbabilitiesFiniteUnphased<4>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 8)
		{
			return getCompressedProbabilitiesFiniteUnphased<8>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else if(nFounders == 16)
		{
			return getCompressedProbabilitiesFiniteUnphased<16>(r, nFunnels, intercrossingGenerations, selfingGenerations);
		}
		else 
		{
			throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
		}
	}
END_RCPP
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
		return getCompressedProbabilitiesInfinite<2>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 4)
	{
		return getCompressedProbabilitiesInfinite<4>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 8)
	{
		return getCompressedProbabilitiesInfinite<8>(r, nFunnels, intercrossingGenerations);
	}
	else if(nFounders == 16)
	{
		return getCompressedProbabilitiesInfinite<16>(r, nFunnels, intercrossingGenerations);
	}
	else 
	{
		throw std::runtime_error("Input nFounders must be one of 2, 4, 8 or 16");
	}
END_RCPP
}
