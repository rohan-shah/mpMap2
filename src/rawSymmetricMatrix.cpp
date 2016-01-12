#include "rawSymmetricMatrix.h"
#include "matrixChunks.h"
SEXP rawSymmetricMatrixSubsetIndices(SEXP object_, SEXP i_, SEXP j_, SEXP drop_)
{
BEGIN_RCPP
	Rcpp::S4 object = object_;
	Rcpp::CharacterVector markers = object.slot("markers");
	R_xlen_t nMarkers = markers.size();
	Rcpp::NumericVector levels = object.slot("levels");
	Rcpp::RawVector data = object.slot("data");
	Rcpp::IntegerVector i = i_;
	Rcpp::IntegerVector j = j_;
	bool drop = Rcpp::as<bool>(drop_);
	if(drop)
	{
		if(i.size() == 1 && j.size() == 1)
		{
			R_xlen_t i = Rcpp::as<int>(i_);
			R_xlen_t j = Rcpp::as<int>(j_);
			if(i > j) std::swap(i, j);
			Rbyte rawValue = data[(j*(j-(R_xlen_t)1))/(R_xlen_t)2 + i-(R_xlen_t)1];
			if(rawValue == 0xff) return Rcpp::wrap(NA_REAL);
			return Rcpp::wrap(levels[rawValue]);
		}
		else if(i.size() == 1)
		{
			Rcpp::NumericVector result(j.size());
			R_xlen_t i = Rcpp::as<int>(i_);
			Rcpp::CharacterVector names(j.size());
			for(R_xlen_t jCounter = 0; jCounter < j.size(); jCounter++)
			{
				R_xlen_t iCopied = i;
				R_xlen_t jValue = j[jCounter];
				names[jCounter] = markers[jValue-(R_xlen_t)1];
				if(iCopied > jValue) std::swap(iCopied, jValue);
				Rbyte rawValue = data[(jValue*(jValue-(R_xlen_t)1))/(R_xlen_t)2 + iCopied-(R_xlen_t)1];
				if(rawValue == 0xff) result[jCounter] = NA_REAL;
				else result[jCounter] = levels[rawValue];
			}
			result.attr("names") = names;
			return result;
		}
		else if(j.size() == 1)
		{
			Rcpp::NumericVector result(i.size());
			R_xlen_t j = Rcpp::as<int>(j_);
			Rcpp::CharacterVector names(i.size());
			for(R_xlen_t iCounter = 0; iCounter < i.size(); iCounter++)
			{
				R_xlen_t jCopied = j;
				R_xlen_t iValue = i[iCounter];
				names[iCounter] = markers[iValue-(R_xlen_t)1];
				if(iValue > jCopied) std::swap(iValue, jCopied);
				Rbyte rawValue = data[(jCopied*(jCopied-(R_xlen_t)1))/(R_xlen_t)2 + iValue-(R_xlen_t)1];
				if(rawValue == 0xff) result[iCounter] = NA_REAL;
				else result[iCounter] = levels[rawValue];
			}
			result.attr("names") = names;
			return result;
		}
	}
	Rcpp::NumericMatrix result((int)i.size(), (int)j.size());
	Rcpp::CharacterVector rownames(i.size()), colnames(j.size());
	for(R_xlen_t iCounter = 0; iCounter < i.size(); iCounter++)
	{
		rownames[iCounter] = markers[i[iCounter]-(R_xlen_t)1];
		for(R_xlen_t jCounter = 0; jCounter < j.size(); jCounter++)
		{
			R_xlen_t iCopied = i[iCounter], jCopied = j[jCounter];
			if(iCopied > jCopied) std::swap(iCopied, jCopied);
			Rbyte rawValue = data[(jCopied*(jCopied-(R_xlen_t)1))/(R_xlen_t)2 + iCopied-(R_xlen_t)1];
			if(rawValue == 0xff) result(iCounter, jCounter) = NA_REAL;
			else result(iCounter, jCounter) = levels[rawValue];
		}
	}
	for(R_xlen_t jCounter = 0; jCounter < j.size(); jCounter++)
	{
		colnames[jCounter] = markers[j[jCounter]-(R_xlen_t)1];
	}
	result.attr("dimnames") = Rcpp::List::create(rownames, colnames);
	return result;
END_RCPP
}
SEXP rawSymmetricMatrixSubsetObject(SEXP object_, SEXP indices_)
{
BEGIN_RCPP
	Rcpp::S4 object = object_;
	Rcpp::RawVector oldData = object.slot("data");
	R_xlen_t oldNMarkers = Rcpp::as<Rcpp::CharacterVector>(object.slot("markers")).size();
	Rcpp::IntegerVector indices = indices_;
	R_xlen_t newNMarkers = indices.size();
	Rcpp::RawVector newData((indices.size() * (indices.size() + (R_xlen_t)1))/(R_xlen_t)2);
	R_xlen_t counter = 0;
	//Column
	for(R_xlen_t j = 0; j < newNMarkers; j++)
	{
		//Row
		for(R_xlen_t i = 0; i <= j; i++)
		{
			R_xlen_t indexJ = indices[j], indexI = indices[i];
			if(indexI > indexJ) std::swap(indexI, indexJ);
			newData(counter) = oldData[(indexJ*(indexJ-(R_xlen_t)1))/(R_xlen_t)2 + indexI - (R_xlen_t)1];
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
	R_xlen_t counter = 0;
	for(; !iterator.isDone(); iterator.next())
	{
		std::pair<int, int> markerPair = iterator.get();
		R_xlen_t markerRow = markerPair.first, markerColumn = markerPair.second;
		destinationData((markerColumn*(markerColumn-(R_xlen_t)1))/(R_xlen_t)2 + (markerRow - (R_xlen_t)1)) = source(counter);
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

	if((indices.size()*(indices.size()+(R_xlen_t)1))/(R_xlen_t)2 != source.size())
	{
		throw std::runtime_error("Mismatch between index length and source object size");
	}
	for(R_xlen_t column = 0; column < indices.size(); column++)
	{
		for(R_xlen_t row = 0; row <= column; row++)
		{
			R_xlen_t rowIndex = indices[row];
			R_xlen_t columnIndex = indices[column];
			if(rowIndex > columnIndex)
			{
				std::swap(rowIndex, columnIndex);
			}
			destinationData((columnIndex*(columnIndex-(R_xlen_t)1))/(R_xlen_t)2+rowIndex-(R_xlen_t)1) = source((column*(column+(R_xlen_t)1))/(R_xlen_t)2 + row);
		}
	}
END_RCPP
}
