#include "arsaRawREntryPoint.h"
#include "arsaRaw.h"
#include <Rcpp.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
SEXP arsaRawREntryPoint(SEXP n_, SEXP rawDist_, SEXP levels_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_, SEXP maxMove_sexp, SEXP effortMultiplier_sexp, SEXP randomStart_sexp)
{
BEGIN_RCPP
	R_xlen_t n;
	try
	{
		n = Rcpp::as<int>(n_);
	}
	catch(...)
	{
		throw std::runtime_error("Input n must be an integer");
	}
	if(n < 1)
	{
		throw std::runtime_error("Input n must be positive");
	}

	Rcpp::RawVector rawDist;
	try
	{
		rawDist = Rcpp::as<Rcpp::RawVector>(rawDist_);
	}
	catch(...)
	{
		throw std::runtime_error("Input dist must be a numeric vector");
	}

	std::vector<double> levels;
	try
	{
		levels = Rcpp::as<std::vector<double> >(levels_);
	}
	catch(...)
	{
		throw std::runtime_error("Input levels must be a numeric vector");
	}

	int nReps;
	try
	{
		nReps = Rcpp::as<int>(nReps_);
	}
	catch(...)
	{
		throw std::runtime_error("Input nReps must be an integer");
	}

	int maxMove;
	try
	{
		maxMove = Rcpp::as<int>(maxMove_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input maxMove must be an integer");
	}
	if(maxMove < 0)
	{
		throw std::runtime_error("Input maxMove must be non-negative");
	}

	bool randomStart;
	try
	{
		randomStart = Rcpp::as<bool>(randomStart_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input randomStart must be a logical");
	}

	double effortMultiplier;
	try
	{
		effortMultiplier = Rcpp::as<double>(effortMultiplier_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input effortMultiplier must be numeric");
	}
	if(effortMultiplier <= 0)
	{
		throw std::runtime_error("Input effortMultiplier must be positive");
	}

	double temperatureMin;
	try
	{
		temperatureMin = Rcpp::as<double>(temperatureMin_);
	}
	catch(...)
	{
		throw std::runtime_error("Input temperatureMin must be a number");
	}
	if(temperatureMin <= 0)
	{
		throw std::runtime_error("Input temperatureMin must be positive");
	}

	double cool;
	try
	{
		cool = Rcpp::as<double>(cool_);
	}
	catch(...)
	{
		throw std::runtime_error("Input cool must be a number");
	}

	//We unpack the rawDist data into a symmetric matrix, for the purposes of running the ordering
	std::vector<Rbyte> distMatrix(n*n);
	for(R_xlen_t i = 0; i < n; i++)
	{
		for(R_xlen_t j = 0; j <= i; j++)
		{
			distMatrix[i * n + j] = distMatrix[j * n + i] = rawDist(i *(i + 1) + j);
		}
	}
	std::vector<int> permutation;
	std::function<bool(long,long)> progressFunction = [](long,long){return false;};
	arsaRaw::arsaRawArgs args(levels, permutation);
	args.n = n;
	args.rawDist = &(distMatrix[0]);
	args.cool = cool;
	args.temperatureMin = temperatureMin;
	args.nReps = nReps;
	args.progressFunction = progressFunction;
	args.randomStart = randomStart;
	args.maxMove = maxMove;
	args.effortMultiplier = effortMultiplier;
#ifdef USE_OPENMP
	if(omp_get_max_threads() > 1)
	{
		arsaRaw::arsaRawParallel(args);
	}
	else
#endif
	{
		arsaRaw::arsaRaw(args);
	}
	return Rcpp::wrap(permutation);
END_RCPP
}
