#include "arsaRaw.h"
#include <Rcpp.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
inline bool descendingComparer(double i, double j)
{
	return i > j;
}
inline void getPairForSwap(R_xlen_t n, R_xlen_t& swap1, R_xlen_t& swap2)
{
	do
	{
		swap1 = (R_xlen_t)(unif_rand()*n);
		swap2 = (R_xlen_t)(unif_rand()*n);
		if(swap1 == n) swap1--;
		if(swap2 == n) swap2--;
	}
	while(swap1 == swap2);
}
inline double deltaFromComponents(std::vector<double>& levels, std::vector<int>& deltaComponents)
{
	double delta = 0;
	for(int i = 0; i < levels.size(); i++)
	{
		delta += deltaComponents[i]*levels[i];
	}
	return delta;
}
inline double computeDelta(std::vector<int>& randomPermutation, R_xlen_t swap1, R_xlen_t swap2, Rbyte* rawDist, std::vector<double>& levels, std::vector<int>& deltaComponents)
{
	R_xlen_t permutationSwap1 = randomPermutation[swap1];
	R_xlen_t permutationSwap2 = randomPermutation[swap2];
	R_xlen_t n = randomPermutation.size();
	std::fill(deltaComponents.begin(), deltaComponents.end(), 0);
	//compute delta
	for(R_xlen_t i = 0; i < n; i++)
	{
		if(i == swap1 || i == swap2) continue;
		R_xlen_t permutationI = randomPermutation[i];
		int count = (int)(abs(i - swap1) - abs(i - swap2));
		deltaComponents[rawDist[permutationSwap2 * n + permutationI]] += count;
		deltaComponents[rawDist[permutationSwap1 * n + permutationI]] -= count;
	}
	//subtract off the case i == swap1.
	//if(permutationSwap2 < permutationSwap1) std::swap(permutationSwap1, permutationSwap2);
	//delta += abs(swap1 - swap2) * dist[(permutationSwap2 * (permutationSwap2+1))/2 + permutationSwap1];
	return deltaFromComponents(levels, deltaComponents);
}
SEXP arsaRaw(SEXP n_, SEXP rawDist_, SEXP levels_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_)
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

	double temperatureMin;
	try
	{
		temperatureMin = Rcpp::as<double>(temperatureMin_);
	}
	catch(...)
	{
		throw std::runtime_error("Input temperatureMin must be a number");
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
	std::function<void(long,long)> progressFunction = [](long,long){};
	arsaRaw((long)n, &(distMatrix[0]), levels, cool, temperatureMin, nReps, permutation, progressFunction);
	return Rcpp::wrap(permutation);
END_RCPP
}
void arsaRaw(long n, Rbyte* rawDist, std::vector<double>& levels, double cool, double temperatureMin, long nReps, std::vector<int>& permutation, std::function<void(long, long)> progressFunction)
{
	permutation.resize(n);
	//We skip the initialisation of D, R1 and R2 from arsa.f, and the computation of asum. 
	//Next the original arsa.f code creates nReps random permutations, and holds them all at once. This doesn't seem necessary, we create them one at a time and discard them
	double zbestAllReps = 0;
	//A copy of the best permutation found
	std::vector<int> bestPermutationThisRep(n);
	//We use this to build the random permutations
	std::vector<int> consecutive(n);
	for(R_xlen_t i = 0; i < n; i++) consecutive[i] = (int)i;
	std::vector<int> deltaComponents(levels.size());
	//We're doing lots of simulation, so we use the old-fashioned approach to dealing with Rs random number generation
	GetRNGstate();

	for(int repCounter = 0; repCounter < nReps; repCounter++)
	{
		//create the random permutation
		for(R_xlen_t i = 0; i < n; i++)
		{
			double rand = unif_rand();
			R_xlen_t index = (R_xlen_t)(rand*(n-i));
			if(index == n-i) index--;
			bestPermutationThisRep[i] = consecutive[index];
			std::swap(consecutive[index], *(consecutive.rbegin()+i));
		}
		//calculate value of z
		double z = 0;
		for(R_xlen_t i = 0; i < n-1; i++)
		{
			R_xlen_t k = bestPermutationThisRep[i];
			for(R_xlen_t j = i+1; j < n; j++)
			{
				R_xlen_t l = bestPermutationThisRep[j];
				z += (j-i) * levels[rawDist[l*n + k]];
			}
		}
		double zbestThisRep = z;
		double temperatureMax = 0;
		//Now try 5000 random swaps
		for(R_xlen_t swapCounter = 0; swapCounter < 5000; swapCounter++)
		{
			R_xlen_t swap1, swap2;
			getPairForSwap(n, swap1, swap2);
			double delta = computeDelta(bestPermutationThisRep, swap1, swap2, rawDist, levels, deltaComponents);
			if(delta < 0)
			{
				if(fabs(delta) > temperatureMax) temperatureMax = fabs(delta);
			}
		}
		double temperature = temperatureMax;
		std::vector<int> currentPermutation = bestPermutationThisRep;
		int nloop = (int)((log(temperatureMin) - log(temperatureMax)) / log(cool));
		long totalSteps = nloop * 100 * n;
		long done = 0;
		long threadZeroCounter = 0;
		//Rcpp::Rcout << "Steps needed: " << nloop << std::endl;
		for(R_xlen_t idk = 0; idk < nloop; idk++)
		{
			//Rcpp::Rcout << "Temp = " << temperature << std::endl;
			for(R_xlen_t k = 0; k < 100*n; k++)
			{
				R_xlen_t swap1, swap2;
				getPairForSwap(n, swap1, swap2);
				//swap
				if(unif_rand() <= 0.5)
				{
					double delta = computeDelta(currentPermutation, swap1, swap2, rawDist, levels, deltaComponents);
					if(delta > -1e-8)
					{
						z += delta;
						std::swap(currentPermutation[swap1], currentPermutation[swap2]);
						if(z > zbestThisRep)
						{
							zbestThisRep = z;
							bestPermutationThisRep = currentPermutation;
						}
					}
					else
					{
						if(unif_rand() <= exp(delta / temperature))
						{
							z += delta;
							std::swap(currentPermutation[swap1], currentPermutation[swap2]);
						}
					}
				}
				//insertion
				else
				{
					//three different patrs of delta
					std::fill(deltaComponents.begin(), deltaComponents.end(), 0);
					int span = (int)abs(swap1 - swap2);
					int span2 = span + 1;
					int permutedSwap1 = currentPermutation[swap1];
					if(swap2 > swap1)
					{
						//compute delta1
						for(R_xlen_t counter1 = swap1+1; counter1 <= swap2; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							for(R_xlen_t counter2 = swap2+1; counter2 < n; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2];
								deltaComponents[rawDist[permutedCounter2*n + permutedCounter1]]++;
							}
							for(R_xlen_t counter2 = 0; counter2 < swap1; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2];
								deltaComponents[rawDist[permutedCounter2*n + permutedCounter1]]--;
							}
						}
						//compute delta2
						for(R_xlen_t counter1 = 0; counter1 < swap1; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] += span;
						}
						for(R_xlen_t counter1 = swap2+1; counter1 < n; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] -= span;
						}
						//compute delta3
						for(R_xlen_t counter1 = swap1+1; counter1 <= swap2; counter1++)
						{
							span2 -= 2;
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] += span2;;
						}
					}
					else
					{
						//compute delta1
						for(R_xlen_t counter1 = swap2; counter1 < swap1; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							for(R_xlen_t counter2 = swap1+1; counter2 < n; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2];
								deltaComponents[rawDist[permutedCounter2*n + permutedCounter1]]--;
							}
							for(R_xlen_t counter2 = 0; counter2 < swap2; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2];
								deltaComponents[rawDist[permutedCounter2*n + permutedCounter1]]++;
							}
						}
						//compute delta2
						for(R_xlen_t counter1 = 0; counter1 < swap2; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] -= span;
						}
						for(R_xlen_t counter1 = swap1+1; counter1 < n; counter1++)
						{
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] += span;
						}
						//compute delta3
						for(R_xlen_t counter1 = swap2; counter1 < swap1; counter1++)
						{
							span2 -= 2;
							R_xlen_t permutedCounter1 = currentPermutation[counter1];
							deltaComponents[rawDist[permutedSwap1*n + permutedCounter1]] -= span2;
						}
					}
					double delta = deltaFromComponents(levels, deltaComponents);
					if(delta > -1e-8 || unif_rand() <= exp(delta / temperature))
					{
						z += delta;
						if(swap2 > swap1)
						{
							for(R_xlen_t i = swap1; i < swap2; i++)
							{
								currentPermutation[i] = currentPermutation[i+1];
							}
							currentPermutation[swap2] = (int)permutedSwap1;
						}
						else
						{
							for(R_xlen_t i = swap1; i > swap2; i--)
							{
								currentPermutation[i] = currentPermutation[i-1];
							}
							currentPermutation[swap2] = (int)permutedSwap1; 
						}
					}
					if(delta > -1e-8 && z > zbestThisRep)
					{
						bestPermutationThisRep = currentPermutation;
						zbestThisRep = z;
					}
				}
#ifdef USE_OPENMY
				#pragma omp atomic
#endif
				{
					done++;
				}
#ifdef USE_OPENMP
				if(omp_get_thread_num() == 0)
#endif
				{
					threadZeroCounter++;
					if(threadZeroCounter % 100 == 0)
					{
						progressFunction(done, totalSteps);
					}
				}

			}
			temperature *= cool;
		}
		if(zbestThisRep > zbestAllReps)
		{
			zbestAllReps = zbestThisRep;
			permutation.swap(bestPermutationThisRep);
		}
	}
	PutRNGstate();
}

