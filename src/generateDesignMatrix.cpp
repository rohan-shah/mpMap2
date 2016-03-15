#include "generateDesignMatrix.h"
SEXP generateDesignMatrix(SEXP n_, SEXP maxOffset_)
{
BEGIN_RCPP
	int n;
	try
	{
		n = Rcpp::as<int>(n_);
	}
	catch(...)
	{
		throw std::runtime_error("Input n must be an integer");
	}

	int maxOffset;
	try
	{
		maxOffset = Rcpp::as<int>(maxOffset_);
	}
	catch(...)
	{
		throw std::runtime_error("Input maxOffset must be an integer");
	}
	if(maxOffset > n)
	{
		throw std::runtime_error("Input maxOffset cannot be larger than n");
	}
	int resultRows = n*maxOffset - maxOffset*(maxOffset - 1)/2;
	Rcpp::IntegerMatrix result(resultRows, n);
	memset(&(result(0, 0)), 0, sizeof(int) * n * resultRows);
	
	for(int i = 0; i < n; i++)
	{
		int offset = 0;
		//j is the section going by rows
		for(int j = 0; j < maxOffset; j++)
		{
			int end = std::min(n - j-1, i) + 1;
			for(int k = std::max(0, i-j); k < end; k++) result(offset + k, i) = 1;
			offset = offset + (n-j);
		}
	}
	return result;
END_RCPP
}
