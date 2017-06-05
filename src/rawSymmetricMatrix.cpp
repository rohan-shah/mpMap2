#include "rawSymmetricMatrix.h"
#include "matrixChunks.h"
SEXP rawSymmetricMatrixUncompress(SEXP object_)
{
BEGIN_RCPP
	Rcpp::S4 object;
	try
	{
		object = object_;
	}
	catch(...)
	{
		throw std::runtime_error("Input object must be an S4 object");
	}

	Rcpp::RawVector data;
	try
	{
		data = Rcpp::as<Rcpp::RawVector>(object.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@data must be a raw vector");
	}

	Rcpp::NumericVector levels;
	try
	{
		levels = Rcpp::as<Rcpp::NumericVector>(object.slot("levels"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@levels must be a numeric vector");
	}
	Rcpp::CharacterVector markers = object.slot("markers");

	int dimension = (int)markers.size();
	Rcpp::NumericMatrix output(dimension, dimension);
	for(int row = 0; row < dimension; row++)
	{
		for(int column = 0; column < dimension; column++)
		{
			R_xlen_t i = row+1;
			R_xlen_t j = column+1;
			if(i > j) std::swap(i, j);
			Rbyte rawValue = data[(j*(j-(R_xlen_t)1))/(R_xlen_t)2 + i-(R_xlen_t)1];
			if(rawValue == 0xff) output(row, column) = NA_REAL;
			else output(row, column) = levels[rawValue];
		}
	}
	return output;
END_RCPP
}
SEXP rawSymmetricMatrixSubsetByMatrix(SEXP object_, SEXP index_)
{
BEGIN_RCPP
	Rcpp::S4 object;
	try
	{
		object = object_;
	}
	catch(...)
	{
		throw std::runtime_error("Input object must be an S4 object");
	}

	Rcpp::RawVector data;
	try
	{
		data = Rcpp::as<Rcpp::RawVector>(object.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@data must be a raw vector");
	}

	Rcpp::NumericVector levels;
	try
	{
		levels = Rcpp::as<Rcpp::NumericVector>(object.slot("levels"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@levels must be a numeric vector");
	}

	Rcpp::IntegerMatrix index;
	try
	{
		index = index_;
	}
	catch(...)
	{
		throw std::runtime_error("Input index must be an integer matrix");
	}
	if(index.ncol() != 2) throw std::runtime_error("Matrix of indices must have two columns");

	int nIndices = index.nrow();
	Rcpp::NumericVector output(nIndices);
	for(int row = 0; row < nIndices; row++)
	{
		R_xlen_t i = index(row, 0);
		R_xlen_t j = index(row, 1);
		if(i > j) std::swap(i, j);
		Rbyte rawValue = data[(j*(j-(R_xlen_t)1))/(R_xlen_t)2 + i-(R_xlen_t)1];
		if(rawValue == 0xff) output(row) = NA_REAL;
		else output(row) = levels[rawValue];
	}
	return output;
END_RCPP
}
SEXP rawSymmetricMatrixSubsetIndices(SEXP object_, SEXP i_, SEXP j_, SEXP drop_)
{
BEGIN_RCPP
	Rcpp::S4 object = object_;
	Rcpp::CharacterVector markers = object.slot("markers");
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
	if(countValuesToEstimate(markerRows, markerColumns) != (unsigned long long)source.size())
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
//Returns the value of any((object@data >= length(object@levels)) & object@data != as.raw(255))
SEXP checkRawSymmetricMatrix(SEXP rawSymmetric_)
{
BEGIN_RCPP
	Rcpp::S4 rawSymmetric = rawSymmetric_;
	Rcpp::NumericVector levels = Rcpp::as<Rcpp::NumericVector>(rawSymmetric.slot("levels"));
	Rcpp::RawVector data = Rcpp::as<Rcpp::RawVector>(rawSymmetric.slot("data"));
	R_xlen_t size = data.size(), levelsSize = levels.size();
	for(R_xlen_t i = 0; i < size; i++)
	{
		if(data[i] >= levelsSize && data[i] != 0xff) return Rcpp::wrap(true);
	}
	return Rcpp::wrap(false);
END_RCPP
}
SEXP rawSymmetricMatrixToDist(SEXP object)
{
BEGIN_RCPP
	Rcpp::S4 rawSymmetric = object;
	Rcpp::NumericVector levels = Rcpp::as<Rcpp::NumericVector>(rawSymmetric.slot("levels"));
	Rcpp::CharacterVector markers = Rcpp::as<Rcpp::CharacterVector>(rawSymmetric.slot("markers"));
	Rcpp::RawVector data = Rcpp::as<Rcpp::RawVector>(rawSymmetric.slot("data"));
	R_xlen_t size = markers.size();

	Rcpp::NumericVector result(size*(size - 1)/2, 0.0);
	int counter = 0;
	for(R_xlen_t row = 0; row < size; row++)
	{
		for(R_xlen_t column = row+1; column < size; column++)
		{
			int byte = data(column * (column + (R_xlen_t)1)/(R_xlen_t)2 + row);
			if(byte == 255) result(counter) = std::numeric_limits<double>::quiet_NaN();
			else result(counter) = levels(byte);
			counter++;
		}
	}
	result.attr("Size") = (int)size;
	result.attr("Labels") = markers;
	result.attr("Diag") = false;
	result.attr("Upper") = false;
	result.attr("class") = "dist";
	return result;
END_RCPP
}
SEXP constructDissimilarityMatrixInternal(unsigned char* data, std::vector<double>& levels, int size, SEXP clusters_, int start, const std::vector<int>& currentPermutation)
{
	Rcpp::IntegerVector clusters = Rcpp::as<Rcpp::IntegerVector>(clusters_);
	int minCluster = *std::min_element(clusters.begin(), clusters.end()), maxCluster = *std::max_element(clusters.begin(), clusters.end());
	if(minCluster != 1)
	{
		throw std::runtime_error("Clusters must have consecutive indices starting at 1");
	}
	std::vector<std::vector<int> > groupIndices(maxCluster);
	for(int i = 0; i < clusters.size(); i++)
	{
		groupIndices[clusters[i]-1].push_back(currentPermutation[i + start]);
	}

	std::vector<int> table(levels.size());

	Rcpp::NumericMatrix result(maxCluster, maxCluster);
	for(int rowCluster = 1; rowCluster <= maxCluster; rowCluster++)
	{
		for(int columnCluster = 1; columnCluster <= rowCluster; columnCluster++)
		{
			const std::vector<int>& columnIndices = groupIndices[columnCluster-1];
			const std::vector<int>& rowIndices = groupIndices[rowCluster-1];
			std::fill(table.begin(), table.end(), 0);
			for(std::vector<int>::const_iterator columnMarker = columnIndices.begin(); columnMarker != columnIndices.end(); columnMarker++)
			{
				for(std::vector<int>::const_iterator rowMarker = rowIndices.begin(); rowMarker != rowIndices.end(); rowMarker++)
				{
					int x = *rowMarker, y = *columnMarker;
					if(x < y) std::swap(x, y);
					int byte = data[x *(x + (R_xlen_t)1)/(R_xlen_t)2 + y];
					if(byte == 255) throw std::runtime_error("Values of NA not allowed");
					table[byte]++;
				}
			}
			double sum = 0;
			for(int i = 0; i < (int)table.size(); i++) sum += table[i] * levels[i];
			result(rowCluster-1, columnCluster-1) = result(columnCluster-1, rowCluster-1) = sum / (columnIndices.size() * rowIndices.size());
		}
	}
	return result;
}
SEXP constructDissimilarityMatrix(SEXP object, SEXP clusters_)
{
BEGIN_RCPP
	Rcpp::S4 rawSymmetric = object;
	Rcpp::NumericVector levels = Rcpp::as<Rcpp::NumericVector>(rawSymmetric.slot("levels"));
	Rcpp::CharacterVector markers = Rcpp::as<Rcpp::CharacterVector>(rawSymmetric.slot("markers"));
	Rcpp::RawVector data = Rcpp::as<Rcpp::RawVector>(rawSymmetric.slot("data"));
	int nMarkers = markers.size();
	std::vector<double> levelsCopied = Rcpp::as<std::vector<double> >(levels);
	
	std::vector<int> permutation(nMarkers);
	for(int i = 0; i < nMarkers; i++) permutation[i] = i;
	return constructDissimilarityMatrixInternal(&(data(0)), levelsCopied, nMarkers, clusters_, 0, permutation);
END_RCPP
}
