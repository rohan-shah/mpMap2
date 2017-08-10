#include "computeAllEpistaticChiSquared.h"
SEXP computeAllEpistaticChiSquared(SEXP probabilities_sexp, SEXP nFounders_sexp, SEXP infiniteSelfing_sexp)
{
BEGIN_RCPP
	Rcpp::S4 probabilities;
	try
	{
		probabilities = Rcpp::as<Rcpp::S4>(probabilities_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input probabilities must be an S4 object");
	}

	Rcpp::NumericMatrix data;
	try
	{
		data = Rcpp::as<Rcpp::NumericMatrix>(probabilities.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Input probabilities@data must be a numeric matrix");
	}

	Rcpp::IntegerMatrix key;
	try
	{
		key = Rcpp::as<Rcpp::IntegerMatrix>(probabilities.slot("key"));
	}
	catch(...)
	{
		throw std::runtime_error("Input probabilities@key must be a numeric matrix");
	}

	bool infiniteSelfing;
	try
	{
		infiniteSelfing = Rcpp::as<bool>(infiniteSelfing_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input infiniteSelfing must be a boolean");
	}

	int nFounders;
	try
	{
		nFounders = Rcpp::as<int>(nFounders_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input nFounders must be a boolean");
	}

	std::vector<int> allAlleles;
	for(int i = 0; i < key.nrow(); i++)
	{
		allAlleles.push_back(key(i, 2));
	}
	std::sort(allAlleles.begin(), allAlleles.end());
	allAlleles.erase(std::unique(allAlleles.begin(), allAlleles.end()), allAlleles.end());

	int nAlleles;
	if(infiniteSelfing) nAlleles = nFounders;
	else nAlleles = (int)allAlleles.size();
	if(data.nrow() % nAlleles != 0) throw std::runtime_error("Input probabilities@data had the wrong number of rows");

	int nMarkers = data.ncol();
	int nLines = data.nrow() / nAlleles;
	Rcpp::NumericMatrix returnValue(nMarkers, nMarkers);
#ifdef USE_OPENMP
	#pragma omp parallel
#endif
	{
		std::vector<double> obs(nAlleles*nAlleles), sums1(nAlleles), sums2(nAlleles);
#ifdef USE_OPENMP
		#pragma omp for schedule(dynamic)
#endif
		for(int marker1 = 0; marker1 < nMarkers; marker1++)
		{
			for(int marker2 = 0; marker2 <= marker1; marker2++)
			{
				std::fill(obs.begin(), obs.end(), 0);
				std::fill(sums1.begin(), sums1.end(), 0);
				std::fill(sums2.begin(), sums2.end(), 0);
				for(int allele1 = 0; allele1 < nAlleles; allele1++)
				{
					for(int lineCounter = 0; lineCounter < nLines; lineCounter++)
					{
						sums1[allele1] += data(lineCounter * nAlleles + allele1, marker1);
						sums2[allele1] += data(lineCounter * nAlleles + allele1, marker2);
					}
					for(int allele2 = 0; allele2 < nAlleles ; allele2++)
					{
						for(int lineCounter = 0; lineCounter < nLines; lineCounter++)
						{
							obs[allele1*nAlleles + allele2] += data(lineCounter * nAlleles + allele1, marker1) * data(lineCounter * nAlleles + allele2, marker2);
						}
					}
				}
				double accumulated = 0;
#ifdef INTERNAL_CHECKS
				double sumExpected = 0, sumObs = 0;
#endif
				for(int allele1 = 0; allele1 < nAlleles; allele1++)
				{
					for(int allele2 = 0; allele2 < nAlleles; allele2++)
					{
						double expected = sums1[allele1] * sums2[allele2] / nLines;
						if(expected != 0)
						{
							double absDifference = fabs(expected - obs[allele1*nAlleles + allele2]);
							accumulated += (absDifference - 0.5) * (absDifference - 0.5) / expected;
						}
#ifdef INTERNAL_CHECKS
						sumExpected += expected;
						sumObs += obs[allele1*nAlleles + allele2];
#endif
					}
				}
#ifdef INTERNAL_CHECKS
				if(fabs(sumExpected - sumObs) > 1e-4) throw std::runtime_error("Internal error");
#endif
				returnValue(marker1, marker2) = returnValue(marker2, marker1) = accumulated;
			}
		}
	}
	return returnValue;
END_RCPP
}
