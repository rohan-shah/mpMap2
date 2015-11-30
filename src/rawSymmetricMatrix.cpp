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
			Rbyte rawValue = data[(j*(j-1))/2 + i-1];
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
				Rbyte rawValue = data[(jValue*(jValue-1))/2 + iCopied-1];
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
				Rbyte rawValue = data[(jCopied*(jCopied-1))/2 + iValue-1];
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
			Rbyte rawValue = data[(jCopied*(jCopied-1))/2 + iCopied-1];
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
	Rcpp::S4 object = object_;
	Rcpp::RawVector oldData = object.slot("data");
	int oldNMarkers = Rcpp::as<Rcpp::CharacterVector>(object.slot("markers")).size();
	Rcpp::IntegerVector indices = indices_;
	int newNMarkers = indices.size();
	Rcpp::RawVector newData((indices.size() * (indices.size() + 1))/2);
	int counter = 0;
	//Column
	for(int j = 0; j < newNMarkers; j++)
	{
		//Row
		for(int i = 0; i <= j; i++)
		{
			newData(counter) = oldData[(indices[j]*(indices[j]-1))/2 + indices[i] - 1];
			counter++;
		}
	}
	return newData;
END_RCPP
}
