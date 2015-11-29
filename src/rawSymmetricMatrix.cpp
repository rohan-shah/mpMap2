#include "rawSymmetricMatrix.h"
SEXP rawSymmetricMatrixSubsetIndices(SEXP object_, SEXP i_, SEXP j_, SEXP drop_)
{
BEGIN_RCPP
	Rcpp::S4 object = object_;
	Rcpp::CharacterVector markers = object.slot("markers");
	int nMarkers = markers.size();
	Rcpp::NumericVector levels = object.slot("levels");
	Rcpp::RawVector data = object.slot("data");
	Rcpp::IntegerVector i = i_;
	Rcpp::IntegerVector j = j_;
	bool drop = Rcpp::as<bool>(drop_);
	if(drop)
	{
		if(i.size() == 1 && j.size() == 1)
		{
			int i = Rcpp::as<int>(i_);
			int j = Rcpp::as<int>(j_);
			if(i > j) std::swap(i, j);
			int inversedY = nMarkers - i;
			Rbyte rawValue = data[nMarkers*(nMarkers+1)/2 - (inversedY+2)*(inversedY+1)/2 + (j - i)];
			if(rawValue == 0xff) return Rcpp::wrap(NA_REAL);
			return Rcpp::wrap(levels[rawValue]);
		}
		else if(i.size() == 1)
		{
			Rcpp::NumericVector result(j.size());
			int i = Rcpp::as<int>(i_);
			for(int jCounter = 0; jCounter < j.size(); jCounter++)
			{
				int iCopied = i;
				int jValue = j[jCounter];
				if(iCopied > jValue) std::swap(iCopied, jValue);
				int inversedY = nMarkers - iCopied;
				Rbyte rawValue = data[nMarkers*(nMarkers+1)/2 - (inversedY+2)*(inversedY+1)/2 + (jValue - i)];
				if(rawValue == 0xff) result[jCounter] = NA_REAL;
				else result[jCounter] = levels[rawValue];
			}
			return result;
		}
		else if(j.size() == 1)
		{
			Rcpp::NumericVector result(i.size());
			int j = Rcpp::as<int>(j_);
			for(int iCounter = 0; iCounter < i.size(); iCounter++)
			{
				int jCopied = j;
				int iValue = i[iCounter];
				if(iValue > jCopied) std::swap(iValue, jCopied);
				int inversedY = nMarkers - iValue;
				Rbyte rawValue = data[nMarkers*(nMarkers+1)/2 - (inversedY+2)*(inversedY+1)/2 + (jCopied - iValue)];
				if(rawValue == 0xff) result[iCounter] = NA_REAL;
				else result[iCounter] = levels[rawValue];
			}
			return result;
		}
	}
	Rcpp::NumericMatrix result(i.size(), j.size());
	for(int iCounter = 0; iCounter < i.size(); iCounter++)
	{
		for(int jCounter = 0; jCounter < j.size(); jCounter++)
		{
			int iCopied = i[iCounter], jCopied = j[jCounter];
			if(iCopied > jCopied) std::swap(iCopied, jCopied);
			int inversedY = nMarkers - iCopied;
			Rbyte rawValue = data[nMarkers*(nMarkers+1)/2 - (inversedY+2)*(inversedY+1)/2 + (jCopied - iCopied)];
			if(rawValue == 0xff) result(iCounter, jCounter) = NA_REAL;
			else result(iCounter, jCounter) = levels[rawValue];
		}
	}
	return result;
END_RCPP
}
SEXP rawSymmetricMatrixSubsetObject(SEXP object_, SEXP indices_)
{
BEGIN_RCPP
	Rcpp::S4 object = object;
	Rcpp::RawVector oldData = object.slot("data");
	int oldNMarkers = Rcpp::as<Rcpp::CharacterVector>(object.slot("markers")).size();
	int newNMarkers = indices.size();
	Rcpp::IntegerVector indices = indices_;
	Rcpp::RawVector newData((indices.size() * (indices.size() + 1)/2);
	int counter = 0;
	//Row
	for(int i = 0; i < newNMarkers; i++)
	{
		//Column
		for(int j = i; j < newNMarkers; j++)
		{
			int inversedY = nMarkers - indices[i];
			newData(counter) = oldData[nMarkers*(nMarkers+1)/2 - (inversedY+2)*(inversedY+1)/2 + (indices[j] - indices[i])];
			counter++;
		}
	}
	return newData;
END_RCPP
}
