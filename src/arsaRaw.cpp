#include "arsaRaw.h"
#include <Rcpp.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace arsaRaw
{
	void arsaRawExported(arsaRawArgs& args)
	{
#ifdef USE_OPENMP
		if(omp_get_max_threads() > 1)
		{
			arsaRawParallel(args);
		}
		else
#endif
		{
			arsaRaw(args);
		}
	}
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
	inline void getPairForMove(R_xlen_t n, R_xlen_t& swap1, R_xlen_t& swap2, int maxMove)
	{
		do
		{
			swap1 = (R_xlen_t)(unif_rand()*n);
			if(maxMove > 0)
			{
				int minSwap2 = std::max((int)swap1 - maxMove, 0);
				int maxSwap2 = std::min((int)swap1 + maxMove, (int)n);
				swap2 = (R_xlen_t)(minSwap2 + unif_rand()*(maxSwap2 - minSwap2));
			}
			else
			{
				swap2 = (R_xlen_t)(unif_rand()*n);
			}
			if(swap1 == n) swap1--;
			if(swap2 == n) swap2--;
		}
		while(swap1 == swap2);
	}

	inline double deltaFromComponents(const std::vector<double>& levels, std::vector<int>& deltaComponents)
	{
		double delta = 0;
		for(int i = 0; i < (int)levels.size(); i++)
		{
			delta += deltaComponents[i]*levels[i];
		}
		return delta;
	}
	inline double computeDelta(const std::vector<int>& randomPermutation, R_xlen_t swap1, R_xlen_t swap2, const Rbyte* rawDist, const std::vector<double>& levels, std::vector<int>& deltaComponents)
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
	inline double computeMoveDelta(std::vector<int>& deltaComponents, int swap1, int swap2, const std::vector<int>& currentPermutation, const Rbyte* rawDist, R_xlen_t n, const std::vector<double>& levels)
	{
		//three different parts of delta
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
		return deltaFromComponents(levels, deltaComponents);
	}
	void arsaRaw(arsaRawArgs& args)
	{
		long n = args.n;
		Rbyte* rawDist = args.rawDist;
		std::vector<double>& levels = args.levels;
		if(std::find(rawDist, rawDist + n*n, (Rbyte)255) != rawDist+n*n)
		{
			throw std::runtime_error("Internal error: Missing values must be imputed before calling arsaRaw");
		}
		double cool = args.cool;
		double temperatureMin = args.temperatureMin;
		if(temperatureMin <= 0)
		{
			throw std::runtime_error("Input temperatureMin must be positive");
		}
		
		long nReps = args.nReps;
		std::vector<int>& permutation = args.permutation;
		std::function<bool(unsigned long, unsigned long)> progressFunction = args.progressFunction;
		bool randomStart = args.randomStart;

		int maxMove = args.maxMove;
		if(maxMove < 0)
		{
			throw std::runtime_error("Input maxMove must be non-negative");
		}

		double effortMultiplier = args.effortMultiplier;
		if(effortMultiplier <= 0)
		{
			throw std::runtime_error("Input effortMultiplier must be positive");
		}

		permutation.resize(n);
		if(n == 1)
		{
			permutation[0] = 0;
			return;
		}
		else if(n < 1)
		{
			throw std::runtime_error("Input n must be positive");
		}
		//We skip the initialisation of D, R1 and R2 from arsa.f, and the computation of asum. 
		//Next the original arsa.f code creates nReps random permutations, and holds them all at once. This doesn't seem necessary, we create them one at a time and discard them
		double zbestAllReps = -std::numeric_limits<double>::infinity();
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
			//create the random permutation, if we decided to use a random initial permutation
			if(randomStart)
			{
				for(R_xlen_t i = 0; i < n; i++)
				{
					double rand = unif_rand();
					R_xlen_t index = (R_xlen_t)(rand*(n-i));
					if(index == n-i) index--;
					bestPermutationThisRep[i] = consecutive[index];
					std::swap(consecutive[index], *(consecutive.rbegin()+i));
				}
			}
			else
			{
				for(R_xlen_t i = 0; i < n; i++)
				{
					bestPermutationThisRep[i] = consecutive[i];
				}
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
			for(R_xlen_t swapCounter = 0; swapCounter < (R_xlen_t)(5000*effortMultiplier); swapCounter++)
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
			long totalSteps = (long)(nloop * 100 * n * effortMultiplier);
			long done = 0;
			long threadZeroCounter = 0;
			//Rcpp::Rcout << "Steps needed: " << nloop << std::endl;
			for(R_xlen_t idk = 0; idk < nloop; idk++)
			{
				//Rcpp::Rcout << "Temp = " << temperature << std::endl;
				for(R_xlen_t k = 0; k < (R_xlen_t)(100*n*effortMultiplier); k++)
				{
					R_xlen_t swap1, swap2;
					//swap
					if(unif_rand() <= 0.5)
					{
						getPairForSwap(n, swap1, swap2);
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
						getPairForMove(n, swap1, swap2, maxMove);
						double delta = computeMoveDelta(deltaComponents, swap1, swap2, currentPermutation, rawDist, n, levels);
						int permutedSwap1 = currentPermutation[swap1];
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
					done++;
					threadZeroCounter++;
					bool cancel = progressFunction(done, totalSteps);
					if(cancel)
					{
						for(int i = 0; i < n; i++) permutation[i] = i;
						PutRNGstate();
						return;
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
#ifdef USE_OPENMP
	//Related to parallel version
	struct change
	{
		bool isMove;
		int swap1, swap2;
		double delta;
	};
	void deltaForChange(change& possibleChange, const std::vector<int>& currentPermutation, const Rbyte* rawDist, const std::vector<double>& levels)
	{
		R_xlen_t swap1 = possibleChange.swap1, swap2 = possibleChange.swap2;
		std::vector<int> deltaComponents(levels.size());
		R_xlen_t n = currentPermutation.size();
		if(possibleChange.isMove)
		{
			possibleChange.delta = computeMoveDelta(deltaComponents, swap1, swap2, currentPermutation, rawDist, n, levels);
		}
		else
		{
			possibleChange.delta = computeDelta(currentPermutation, swap1, swap2, rawDist, levels, deltaComponents);
		}
	}
	void makeChange(change& possibleChange, std::vector<int>& currentPermutation, const Rbyte* rawDist, const std::vector<double>& levels, double z, double zbestThisRep, std::vector<int>& bestPermutationThisRep, double temperature)
	{
		R_xlen_t swap1 = possibleChange.swap1, swap2 = possibleChange.swap2;
		double delta = possibleChange.delta;
		if(possibleChange.isMove)
		{
			int permutedSwap1 = currentPermutation[swap1];
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
		else
		{
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
	}
	void arsaRawParallel(arsaRawArgs& args)
	{
		long n = args.n;
		Rbyte* rawDist = args.rawDist;
		std::vector<double>& levels = args.levels;
		if(std::find(rawDist, rawDist + n*n, (Rbyte)255) != rawDist+n*n)
		{
			throw std::runtime_error("Internal error: Missing values must be imputed before calling arsaRawParallel");
		}
		double cool = args.cool;
		double temperatureMin = args.temperatureMin;
		if(temperatureMin <= 0)
		{
			throw std::runtime_error("Input temperatureMin must be positive");
		}
		
		long nReps = args.nReps;
		std::vector<int>& permutation = args.permutation;
		std::function<bool(unsigned long, unsigned long)> progressFunction = args.progressFunction;
		bool randomStart = args.randomStart;

		int maxMove = args.maxMove;
		if(maxMove < 0)
		{
			throw std::runtime_error("Input maxMove must be non-negative");
		}

		double effortMultiplier = args.effortMultiplier;
		if(effortMultiplier <= 0)
		{
			throw std::runtime_error("Input effortMultiplier must be positive");
		}

		permutation.resize(n);
		if(n == 1)
		{
			permutation[0] = 0;
			return;
		}
		else if(n < 1)
		{
			throw std::runtime_error("Input n must be positive");
		}

		//We skip the initialisation of D, R1 and R2 from arsa.f, and the computation of asum. 
		//Next the original arsa.f code creates nReps random permutations, and holds them all at once. This doesn't seem necessary, we create them one at a time and discard them
		double zbestAllReps = -std::numeric_limits<double>::infinity();
		//A copy of the best permutation found
		std::vector<int> bestPermutationThisRep(n);
		//We use this to build the random permutations
		std::vector<int> consecutive(n);
		for(R_xlen_t i = 0; i < n; i++) consecutive[i] = (int)i;
		std::vector<int> deltaComponents(levels.size());
		//We're doing lots of simulation, so we use the old-fashioned approach to dealing with Rs random number generation
		GetRNGstate();

		std::vector<change> stackOfChanges;
		std::vector<bool> dirty(n, false);
		for(int repCounter = 0; repCounter < nReps; repCounter++)
		{
			//create the random permutation, if we decided to use a random initial permutation
			if(randomStart)
			{
				for(R_xlen_t i = 0; i < n; i++)
				{
					double rand = unif_rand();
					R_xlen_t index = (R_xlen_t)(rand*(n-i));
					if(index == n-i) index--;
					bestPermutationThisRep[i] = consecutive[index];
					std::swap(consecutive[index], *(consecutive.rbegin()+i));
				}
			}
			else
			{
				for(R_xlen_t i = 0; i < n; i++)
				{
					bestPermutationThisRep[i] = consecutive[i];
				}
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
			for(R_xlen_t swapCounter = 0; swapCounter < (R_xlen_t)(5000*effortMultiplier); swapCounter++)
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
			long totalSteps = (long)(nloop * 100 * n * effortMultiplier);
			long done = 0;
			//Rcpp::Rcout << "Steps needed: " << nloop << std::endl;
			for(R_xlen_t idk = 0; idk < nloop; idk++)
			{
				//Rcpp::Rcout << "Temp = " << temperature << std::endl;
				for(R_xlen_t k = 0; k < (R_xlen_t)(100*n*effortMultiplier); k++)
				{
					R_xlen_t swap1, swap2;
					//swap
					if(unif_rand() <= 0.5)
					{
						getPairForSwap(n, swap1, swap2);
						change newChange;
						newChange.isMove = false;
						newChange.swap1 = swap1; newChange.swap2 = swap2;

						if(dirty[swap1] || dirty[swap2])
						{
							#pragma omp parallel for
							for(std::vector<change>::iterator i = stackOfChanges.begin(); i < stackOfChanges.end(); i++)
							{
								deltaForChange(*i, currentPermutation, rawDist, levels);
							}
							for(std::vector<change>::iterator i = stackOfChanges.begin(); i != stackOfChanges.end(); i++)
							{
								makeChange(*i, currentPermutation, rawDist, levels, z, zbestThisRep, bestPermutationThisRep, temperature);
							}
							done += stackOfChanges.size();
							bool cancel = progressFunction(done, totalSteps);
							if(cancel)
							{
								for(int i = 0; i < n; i++) permutation[i] = i;
								PutRNGstate();
								return;
							}
							stackOfChanges.clear();
							std::fill(dirty.begin(), dirty.end(), false);
						}
						else dirty[swap1] = dirty[swap2] = true;
						stackOfChanges.push_back(newChange);
					}
					//insertion
					else
					{
						getPairForMove(n, swap1, swap2, maxMove);
						bool canDefer = true;
						for(R_xlen_t i = std::min(swap1, swap2); i != std::max(swap1, swap2)+1; i++) canDefer &= !dirty[i];
						change newChange;
						newChange.isMove = true;
						newChange.swap1 = swap1; 
						newChange.swap2 = swap2;
						if(canDefer)
						{
							std::fill(dirty.begin() + std::min(swap1, swap2), dirty.begin() + std::max(swap1, swap2)+1, true);
						}
						else
						{
							#pragma omp parallel for
							for(std::vector<change>::iterator i = stackOfChanges.begin(); i < stackOfChanges.end(); i++)
							{
								deltaForChange(*i, currentPermutation, rawDist, levels);
							}
							for(std::vector<change>::iterator i = stackOfChanges.begin(); i != stackOfChanges.end(); i++)
							{
								makeChange(*i, currentPermutation, rawDist, levels, z, zbestThisRep, bestPermutationThisRep, temperature);
							}

							done += stackOfChanges.size();
							bool cancel = progressFunction(done, totalSteps);
							if(cancel)
							{
								for(int i = 0; i < n; i++) permutation[i] = i;
								PutRNGstate();
								return;
							}
							stackOfChanges.clear();
							std::fill(dirty.begin(), dirty.end(), false);
						}
						stackOfChanges.push_back(newChange);
					}
				}
				#pragma omp parallel for
				for(std::vector<change>::iterator i = stackOfChanges.begin(); i < stackOfChanges.end(); i++)
				{
					deltaForChange(*i, currentPermutation, rawDist, levels);
				}
				for(std::vector<change>::iterator i = stackOfChanges.begin(); i != stackOfChanges.end(); i++)
				{
					makeChange(*i, currentPermutation, rawDist, levels, z, zbestThisRep, bestPermutationThisRep, temperature);
				}

				done += stackOfChanges.size();
				bool cancel = progressFunction(done, totalSteps);
				if(cancel)
				{
					for(int i = 0; i < n; i++) permutation[i] = i;
					PutRNGstate();
					return;
				}
				stackOfChanges.clear();
				std::fill(dirty.begin(), dirty.end(), false);
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
#endif
}
