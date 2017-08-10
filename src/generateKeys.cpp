#include "generateKeys.h"
void generateKeys(Rcpp::IntegerMatrix& key, Rcpp::IntegerMatrix& outputKey, int nFounders, bool infiniteSelfing)
{
	key = Rcpp::IntegerMatrix(nFounders, nFounders);
	for(int i = 0; i < nFounders; i++)
	{
		key(i, i) = i + 1;
	}
	int counter = nFounders+1;
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = i+1; j < nFounders; j++)
		{
			key(j, i) = key(i, j) = counter;
			counter++;
		}
	}
	//We also want a version closer to the hetData format
	if(infiniteSelfing)
	{
		outputKey = Rcpp::IntegerMatrix(nFounders, 3);
		for(int i = 0; i < nFounders; i++)
		{
			outputKey(i, 0) = i+1;
			outputKey(i, 1) = i+1;
			outputKey(i, 2) = key(i, i);
		}
	}
	else
	{
		outputKey = Rcpp::IntegerMatrix(nFounders*nFounders, 3);
		int counter = 0;
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				outputKey(counter, 0) = i+1;
				outputKey(counter, 1) = j+1;
				outputKey(counter, 2) = key(i, j);
				counter++;
			}
		}
	}
}
