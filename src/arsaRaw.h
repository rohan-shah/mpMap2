#ifndef MPMAP2_ARSA_RAW_HEADER_GUARD
#define MPMAP2_ARSA_RAW_HEADER_GUARD
#include "Rcpp.h"
#include <functional>
namespace arsaRaw
{
	struct arsaRawArgs
	{
	public:
		arsaRawArgs(std::vector<double>& levels, std::vector<int>& permutation)
			:n(-1), rawDist(NULL), cool(0.5), temperatureMin(0.1), nReps(1), randomStart(true), maxMove(0), effortMultiplier(1), levels(levels), permutation(permutation)
		{}
		long n;
		Rbyte* rawDist;
		double cool;
		double temperatureMin;
		long nReps;
		std::function<bool(unsigned long,unsigned long)> progressFunction;
		bool randomStart;
		int maxMove;
		double effortMultiplier;
		std::vector<double>& levels;
		std::vector<int>& permutation;
	};
	void arsaRaw(arsaRawArgs& args);
	void arsaRawExported(arsaRawArgs& args);
#ifdef _OPENMP
	void arsaRawParallel(arsaRawArgs& args);
#endif
}
#endif
