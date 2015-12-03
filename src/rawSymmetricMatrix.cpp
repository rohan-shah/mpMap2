#include "rawSymmetricMatrix.h"
#include "matrixChunks.h"
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
			int indexJ = indices[j], indexI = indices[i];
			if(indexI > indexJ) std::swap(indexI, indexJ);
			newData(counter) = oldData[(indexJ*(indexJ-1))/2 + indexI - 1];
			counter++;
		}
	}
	return newData;
END_RCPP
}
SEXP assignRawSymmetricMatrixFromEstimateRF(SEXP destination_, SEXP rowIndices_, SEXP columnIndices_, SEXP source_)
{
BEGIN_RCPP
	Rcpp::S4 destination = destination_;
	Rcpp::RawVector source = source_;
	Rcpp::RawVector destinationData = destination.slot("data");
	Rcpp::IntegerVector rowIndices = rowIndices_;
	Rcpp::IntegerVector columnIndices = columnIndices_;

	if(&(source(0)) == &(destinationData(0)))
	{
		throw std::runtime_error("Source and destination cannot be the same in assignRawSymmetricMatrixDiagonal");
	}

	std::vector<int> markerRows, markerColumns;
	markerRows = Rcpp::as<std::vector<int> >(rowIndices);
	markerColumns = Rcpp::as<std::vector<int> >(columnIndices);
	if(countValuesToEstimate(markerRows, markerColumns) != source.size())
	{
		throw std::runtime_error("Mismatch between index length and source object size");
	}

	triangularIterator iterator(markerRows, markerColumns);
	int counter = 0;
	for(; !iterator.isDone(); iterator.next())
	{
		std::pair<int, int> markerPair = iterator.get();
		int markerRow = markerPair.first, markerColumn = markerPair.second;
		destinationData((markerColumn*(markerColumn-1))/2 + (markerRow - 1)) = source(counter);
		counter++;
	}
	return R_NilValue;
END_RCPP
}
SEXP assignRawSymmetricMatrixDiagonal(SEXP destination_, SEXP indices_, SEXP source_)
{
BEGIN_RCPP
	Rcpp::S4 destination = destination_;
	Rcpp::RawVector source = source_;
	Rcpp::RawVector destinationData = destination.slot("data");
	Rcpp::IntegerVector indices = indices_;

	if(&(source(0)) == &(destinationData(0)))
	{
		throw std::runtime_error("Source and destination cannot be the same in assignRawSymmetricMatrixDiagonal");
	}

	if((indices.size()*(indices.size()+1))/2 != source.size())
	{
		throw std::runtime_error("Mismatch between index length and source object size");
	}
	for(int column = 0; column < indices.size(); column++)
	{
		for(int row = 0; row <= column; row++)
		{
			int rowIndex = indices[row];
			int columnIndex = indices[column];
			if(rowIndex > columnIndex)
			{
				std::swap(rowIndex, columnIndex);
			}
			destinationData((columnIndex*(columnIndex-1))/2+rowIndex-1) = source((column*(column+1))/2 + row);
		}
	}
END_RCPP
}
