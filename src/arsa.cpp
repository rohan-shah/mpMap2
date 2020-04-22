#include "arsa.h"
#include <Rcpp.h>
inline bool descendingComparer(double i, double j)
{
	return i > j;
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
inline double computeDelta(std::vector<int>& randomPermutation, R_xlen_t swap1, R_xlen_t swap2, double* dist)
{
	double delta = 0;
	R_xlen_t permutationSwap1 = randomPermutation[swap1];
	R_xlen_t permutationSwap2 = randomPermutation[swap2];
	R_xlen_t n = randomPermutation.size();
	//compute delta
	for(R_xlen_t i = 0; i < n; i++)
	{
		if(i == swap1 || i == swap2) continue;
		R_xlen_t permutationI_1 = randomPermutation[i], copiedPermutationSwap1 = permutationSwap1, copiedPermutationSwap2 = permutationSwap2;
		R_xlen_t permutationI_2 = permutationI_1;
		if(copiedPermutationSwap1 < permutationI_1) std::swap(copiedPermutationSwap1, permutationI_1);
		if(copiedPermutationSwap2 < permutationI_2) std::swap(copiedPermutationSwap2, permutationI_2);
		delta += (std::abs(i - swap1) - std::abs(i - swap2)) * (dist[(copiedPermutationSwap2 * (copiedPermutationSwap2+1))/2 + permutationI_2] - dist[(copiedPermutationSwap1 *(copiedPermutationSwap1 + 1))/2 + permutationI_1]);
	}
	//subtract off the case i == swap1.
	//if(permutationSwap2 < permutationSwap1) std::swap(permutationSwap1, permutationSwap2);
	//delta += std::abs(swap1 - swap2) * dist[(permutationSwap2 * (permutationSwap2+1))/2 + permutationSwap1];
	return delta;
}
SEXP arsaExportedR(SEXP n_, SEXP dist_, SEXP cool_, SEXP temperatureMin_, SEXP nReps_, SEXP maxMove_sexp, SEXP effortMultiplier_sexp, SEXP randomStart_sexp)
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
	if(n == 1)
	{
		Rcpp::IntegerVector retVal(1);
		retVal[0] = 0;
		return retVal;
	}
	else if(n < 1)
	{
		throw std::runtime_error("Input n must be positive. Is this a zero dimension matrix?");
	}

	Rcpp::NumericVector dist;
	try
	{
		dist = Rcpp::as<Rcpp::NumericVector>(dist_);
	}
	catch(...)
	{
		throw std::runtime_error("Input dist must be a numeric vector");
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
	if(cool <= 0)
	{
		throw std::runtime_error("Input cool must be positive");
	}

	arsaArgs args;
	args.n = n;
	args.dist = &(dist(0));
	args.nReps = nReps;
	args.temperatureMin = temperatureMin;
	args.cool = cool;
	args.randomStart = randomStart;
	args.maxMove = maxMove;
	args.effortMultiplier = effortMultiplier;
	std::function<void(unsigned long,unsigned long)> noProgress = [](unsigned long,unsigned long){};
	args.progressFunction = noProgress;
	arsa(args);
	return Rcpp::wrap(args.bestPermutationAllReps);
END_RCPP
}
void arsa(arsaArgs& args)
{
	R_xlen_t n = args.n;
	if(n <= 0)
	{
		throw std::runtime_error("Input n must be positive. Is this a zero dimension matrix?");
	}
	double* dist = args.dist;
	int nReps = args.nReps;
	if(nReps <= 0)
	{
		throw std::runtime_error("Input nReps must be positive");
	}
	double temperatureMin = args.temperatureMin;
	if(temperatureMin <= 0)
	{
		throw std::runtime_error("Input temperatureMin must be positive");
	}
	double cool = args.cool;
	
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
	bool randomStart = args.randomStart;
	std::function<void(unsigned long,unsigned long)> progressFunction = args.progressFunction;
	//We skip the initialisation of D, R1 and R2 from arsa.f, and the computation of asum. 
	//Next the original arsa.f code creates nReps random permutations, and holds them all at once. This doesn't seem necessary, we create them one at a time and discard them
	double zbestAllReps = 0;
	//A copy of the best permutation found
	std::vector<int> bestPermutationThisRep(n);
	std::vector<int>& bestPermutationAllReps = args.bestPermutationAllReps;
	bestPermutationAllReps.resize(n);

	//We use this to build the random permutations
	std::vector<int> consecutive(n);
	for(R_xlen_t i = 0; i < n; i++) consecutive[i] = (int)i;
	//We're doing lots of simulation, so we use the old-fashioned approach to dealing with Rs random number generation
	GetRNGstate();
	bool diagnostics = false;

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
				R_xlen_t kCopied = k;
				if(l < k) std::swap(l, kCopied);
				z += (j-i) * dist[(l*(l+1))/2 + kCopied];
			}
		}
		double zbestThisRep = z;
		double temperatureMax = 0;
		//Now try 5000 random swaps
		for(R_xlen_t swapCounter = 0; swapCounter < (R_xlen_t)(5000*effortMultiplier); swapCounter++)
		{
			R_xlen_t swap1, swap2;
			getPairForSwap(n, swap1, swap2);
			double delta = computeDelta(bestPermutationThisRep, swap1, swap2, dist);
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
		if(diagnostics) Rcpp::Rcout << "Steps needed: " << nloop << std::endl;
		for(R_xlen_t idk = 0; idk < nloop; idk++)
		{
			if(diagnostics) Rcpp::Rcout << "Temp = " << temperature << std::endl;
			for(R_xlen_t k = 0; k < (R_xlen_t)(100*n*effortMultiplier); k++)
			{
				R_xlen_t swap1, swap2;
				//swap
				if(unif_rand() <= 0.5)
				{
					getPairForSwap(n, swap1, swap2);
					double delta = computeDelta(currentPermutation, swap1, swap2, dist);
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
					//three different patrs of delta
					double delta1 = 0, delta2 = 0, delta3 = 0;
					R_xlen_t span = std::abs(swap1 - swap2);
					R_xlen_t span2 = span + 1;
					R_xlen_t permutedSwap1 = currentPermutation[swap1];
					if(swap2 > swap1)
					{
						//compute delta1
						for(R_xlen_t counter1 = swap1+1; counter1 <= swap2; counter1++)
						{
							for(R_xlen_t counter2 = swap2+1; counter2 < n; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2], permutedCounter1 = currentPermutation[counter1];
								if(permutedCounter1 > permutedCounter2) std::swap(permutedCounter1, permutedCounter2);
								//permutedCounter1 = row, permutedCounter2 = column
								delta1 += dist[(permutedCounter2*(permutedCounter2+1))/2 + permutedCounter1];
							}
							for(R_xlen_t counter2 = 0; counter2 < swap1; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2], permutedCounter1 = currentPermutation[counter1];
								if(permutedCounter1 > permutedCounter2) std::swap(permutedCounter1, permutedCounter2);
								//permutedCounter1 = row, permutedCounter2 = column
								delta1 -= dist[(permutedCounter2*(permutedCounter2+1))/2 + permutedCounter1];
							}
						}
						//compute delta2
						for(R_xlen_t counter1 = 0; counter1 < swap1; counter1++)
						{
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(permutedCounter1 > copiedPermutedSwap1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//copiedPermutedSwap1 = column,  permutedCounter1 = row
							delta2 += dist[(copiedPermutedSwap1*(copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
						for(R_xlen_t counter1 = swap2+1; counter1 < n; counter1++)
						{
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(permutedCounter1 > copiedPermutedSwap1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//copiedPermutedSwap1 = column,  permutedCounter1 = row
							delta2 -= dist[(copiedPermutedSwap1*(copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
						//compute delta3
						for(R_xlen_t counter1 = swap1+1; counter1 <= swap2; counter1++)
						{
							span2 -= 2;
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(copiedPermutedSwap1 < permutedCounter1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//permutedCounter1 = row, copiedPermutedSwap1 = column
							delta3 += span2 * dist[(copiedPermutedSwap1 * (copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
					}
					else
					{
						//compute delta1
						for(R_xlen_t counter1 = swap2; counter1 < swap1; counter1++)
						{
							for(R_xlen_t counter2 = swap1+1; counter2 < n; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2], permutedCounter1 = currentPermutation[counter1];
								if(permutedCounter1 > permutedCounter2) std::swap(permutedCounter1, permutedCounter2);
								//permutedCounter1 = row, permutedCounter2 = column
								delta1 -= dist[(permutedCounter2*(permutedCounter2+1))/2 + permutedCounter1];
							}
							for(R_xlen_t counter2 = 0; counter2 < swap2; counter2++)
							{
								R_xlen_t permutedCounter2 = currentPermutation[counter2], permutedCounter1 = currentPermutation[counter1];
								if(permutedCounter1 > permutedCounter2) std::swap(permutedCounter1, permutedCounter2);
								//permutedCounter1 = row, permutedCounter2 = column
								delta1 += dist[(permutedCounter2*(permutedCounter2+1))/2 + permutedCounter1];
							}
						}
						//compute delta2
						for(R_xlen_t counter1 = 0; counter1 < swap2; counter1++)
						{
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(permutedCounter1 > copiedPermutedSwap1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//copiedPermutedSwap1 = column,  permutedCounter1 = row
							delta2 -= dist[(copiedPermutedSwap1*(copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
						for(R_xlen_t counter1 = swap1+1; counter1 < n; counter1++)
						{
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(permutedCounter1 > copiedPermutedSwap1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//copiedPermutedSwap1 = column,  permutedCounter1 = row
							delta2 += dist[(copiedPermutedSwap1*(copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
						//compute delta3
						for(R_xlen_t counter1 = swap2; counter1 < swap1; counter1++)
						{
							span2 -= 2;
							R_xlen_t copiedPermutedSwap1 = permutedSwap1, permutedCounter1 = currentPermutation[counter1];
							if(copiedPermutedSwap1 < permutedCounter1) std::swap(copiedPermutedSwap1, permutedCounter1);
							//permutedCounter1 = row, copiedPermutedSwap1 = column
							delta3 -= span2 * dist[(copiedPermutedSwap1 * (copiedPermutedSwap1 + 1))/2 + permutedCounter1];
						}
					}
					double delta = delta1 + span * delta2 + delta3;
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
				if(threadZeroCounter % 100 == 0)
				{
					progressFunction(done, totalSteps);
				}
			}
			temperature *= cool;
		}
		if(zbestThisRep > zbestAllReps)
		{
			zbestAllReps = zbestThisRep;
			bestPermutationAllReps.swap(bestPermutationThisRep);
		}
	}
	PutRNGstate();
}
