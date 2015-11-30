#include "estimateRF.h"
#include <math.h>
#include <limits>
#include "estimateRFSpecificDesign.h"
#include <stdexcept>
#include "matrixChunks.h"
SEXP estimateRF(SEXP object_, SEXP recombinationFractions_, SEXP markerRows_, SEXP markerColumns_, SEXP lineWeights_, SEXP keepLod_, SEXP keepLkhd_)
{
	BEGIN_RCPP
		Rcpp::NumericVector recombinationFractions;
		try
		{
			recombinationFractions = recombinationFractions_;
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input recombinationFractions must be a numeric vector");
		}
		int nRecombLevels = recombinationFractions.size();
		Rcpp::NumericVector::iterator halfIterator = std::find(recombinationFractions.begin(), recombinationFractions.end(), 0.5);
		if(halfIterator == recombinationFractions.end()) throw std::runtime_error("Input recombinationFractions did not contain the value 0.5");
		int halfIndex = (int)std::distance(recombinationFractions.begin(), halfIterator);

		//Confirm that recombination fractions are in increasing order
		for(int i = 0; i < recombinationFractions.size()-1; i++)
		{
			if(recombinationFractions[i] >= recombinationFractions[i+1]) throw std::runtime_error("Input recombinationFractions must be a vector of increasing values");
		}

		std::vector<double> recombinationFractionsDouble = Rcpp::as<std::vector<double> >(recombinationFractions);
		Rcpp::S4 object;
		try
		{
			object = object_;
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input object must be an S4 object");
		}
		std::vector<int> markerRows;
		try
		{
			markerRows = Rcpp::as<std::vector<int> >(markerRows_);
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input markerRows must be an integer vector");
		}
		std::sort(markerRows.begin(), markerRows.end());
		for(std::vector<int>::iterator markerRow = markerRows.begin(); markerRow != markerRows.end(); markerRow++)
		{
			(*markerRow)--;
		}

		std::vector<int> markerColumns;
		try
		{
			markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input markerColumns must be an integer vector");
		}
		std::sort(markerColumns.begin(), markerColumns.end());
		for(std::vector<int>::iterator markerColumn = markerColumns.begin(); markerColumn != markerColumns.end(); markerColumn++)
		{
			(*markerColumn)--;
		}
		Rcpp::List lineWeights;
		try
		{
			lineWeights = lineWeights_;
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input lineWeights must be a list");
		}
	
		Rcpp::List geneticData;
		try
		{
			geneticData = object.slot("geneticData");
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input object must have a slot named \"geneticData\" which must be a list");
		}
		int nDesigns = geneticData.length();
		if(lineWeights.size() != nDesigns) throw std::runtime_error("Inut lineWeights had the wrong number of entries");
		try
		{
			for(int i = 0; i < nDesigns; i++)
			{
				Rcpp::NumericVector currentDesignLineWeights = lineWeights(i);
			}
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input lineWeights must be a list of numeric vectors");
		}
		bool keepLod, keepLkhd;
		try
		{
			keepLod = Rcpp::as<bool>(keepLod_);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input keepLod must be a boolean");
		}
		try
		{
			keepLkhd = Rcpp::as<bool>(keepLkhd_);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Input keepLkhd must be a boolean");
		}
		if(nDesigns <= 0) throw std::runtime_error("There must be at least one design");
		if(markerRows.size() == 0) throw std::runtime_error("Input markerRows must have at least one entry");
		if(markerColumns.size() == 0) throw std::runtime_error("Input markerColumns must have at least one entry");

		int markerRowMin = *std::min_element(markerRows.begin(), markerRows.end());
		int markerColumnMin = *std::min_element(markerColumns.begin(), markerColumns.end());
		if(markerRowMin < 0) throw std::runtime_error("Invalid values for input markerRows");
		if(markerColumnMin < 0) throw std::runtime_error("Invalid value for input markerColumns");

		//If the input values of markerRows and markerColumns give a region that's completely in the lower triangular region, then throw an error
		unsigned long long nValuesToEstimate = countValuesToEstimate(markerRows, markerColumns);
		if(nValuesToEstimate == 0)
		{
			throw std::runtime_error("Input values of markerRows and markerColumns give a region that is contained in the lower triangular part of the matrix");
		}

		//Warn if we're going to allocate over 4gb
		if(nValuesToEstimate > 4000000000ULL)
		{
			Rcpp::Rcout << "Allocating results matrix of " << (nValuesToEstimate * (std::size_t)nRecombLevels) << " bytes" << std::endl;
		}
		//This is not an Rcpp::NumericVector because it can quite easily overflow the size of such a vector (signed int)
		std::vector<double> result(nValuesToEstimate * nRecombLevels, 0);

		//Last bit of validation
		for(int i = 0; i < nDesigns; i++)
		{
			Rcpp::S4 currentGeneticData = geneticData(i);
			Rcpp::IntegerMatrix finals = currentGeneticData.slot("finals");
			std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
			if((int)lineWeightsThisDesign.size() != finals.nrow())
			{
				throw std::runtime_error("An entry of input lineWeights had the wrong length");
			}
		}
		//Construct vector of argument Objects
		std::vector<estimateRFSpecificDesignArgs> argumentObjects;
		for(int i = 0; i < nDesigns; i++)
		{
			Rcpp::S4 currentGeneticData = geneticData(i);
			std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
			std::string error;
			estimateRFSpecificDesignArgs args(recombinationFractionsDouble, markerRows, markerColumns);
			try
			{
				args.founders = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("founders"));
			}
			catch(Rcpp::not_compatible&)
			{
				std::stringstream ss; 
				ss << "Founders slot of design " << i << " was not an integer matrix";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.finals = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("finals"));
			}
			catch(Rcpp::not_compatible&)
			{
				std::stringstream ss; 
				ss << "Finals slot of design " << i << " was not an integer matrix";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.pedigree = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("pedigree"));
			}
			catch(Rcpp::not_compatible&)
			{
				std::stringstream ss; 
				ss << "Pedigree slot of design " << i << " was not an S4 object";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.hetData = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("hetData"));
			}
			catch(Rcpp::not_compatible&)
			{
				std::stringstream ss; 
				ss << "hetData slot of design " << i << " was not an S4 object";
				throw std::runtime_error(ss.str().c_str());
			}
			args.result = &(result[0]);
			args.lineWeights.swap(lineWeightsThisDesign);
			argumentObjects.emplace_back(args);
		}
		//Estimate required memory usage
		unsigned long long lookupBytes = 0;
		for(int i = 0; i < nDesigns; i++)
		{
			unsigned long long currentLookupBytes = estimateLookup(argumentObjects[i]);
			if(currentLookupBytes == 0) throw std::runtime_error(argumentObjects[i].error.c_str());
			lookupBytes = std::max(lookupBytes, currentLookupBytes);
		}
		if(lookupBytes > 1000000000)
		{
			Rcpp::Rcout << "Maximum lookup table size of " << lookupBytes << " bytes" << std::endl;
		}
		//Now the actual computation
		for(int i = 0; i < nDesigns; i++)
		{
			bool successful = estimateRFSpecificDesign(argumentObjects[i]);
			if(!successful) throw std::runtime_error(argumentObjects[i].error);
		}
		//now for some post-processing to get out the MLE, lod (maybe) and lkhd (maybe)
		Rcpp::NumericVector lod, lkhd;
		Rcpp::RawVector theta(nValuesToEstimate);
		if(keepLod) lod = Rcpp::NumericVector(nValuesToEstimate);
		if(keepLkhd) lkhd = Rcpp::NumericVector(nValuesToEstimate);
	
		double* resultPtr = &(result[0]);
#ifdef USE_OPENMP
		#pragma omp parallel for schedule(static, 1)
#endif
		for(int counter = 0; counter < nValuesToEstimate; counter++)
		{
			double* start = resultPtr + counter * nRecombLevels;
			double* end = resultPtr + (counter + 1) * nRecombLevels;
			double* maxPtr = std::max_element(start, end), *minPtr = std::min_element(start, end);
			double max = *maxPtr, min = *minPtr;
			int currentTheta;
			double currentLod;
			//This is the case where no data was available, across any of the experiments. This is precise, no numerical error involved
			if(max == 0 && min == 0)
			{
				currentLod = std::numeric_limits<double>::quiet_NaN();
				currentTheta = 0xff;
			}
			else
			{
				currentTheta = (int)(maxPtr - start);
				currentLod = max - resultPtr[counter * (unsigned long long)nRecombLevels + halfIndex];
			}
			theta(counter) = currentTheta;
			if(keepLkhd) lkhd(counter) = max;
			if(keepLod) lod(counter) = currentLod;
		}
		Rcpp::RObject lodRet, lkhdRet;
		
		if(keepLod) lodRet = lod;
		else lodRet = R_NilValue;

		if(keepLkhd) lkhdRet = lkhd;
		else lkhdRet = R_NilValue;

		return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lod") = lodRet, Rcpp::Named("lkhd") = lkhdRet, Rcpp::Named("r") = recombinationFractions);
	END_RCPP
}
