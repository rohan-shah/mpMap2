#include "estimateRF.h"
#include <math.h>
#include <limits>
#include "estimateRFSpecificDesign.h"
#include <stdexcept>
#include "matrixChunks.h"
SEXP estimateRF(SEXP object_, SEXP recombinationFractions_, SEXP markerRows_, SEXP markerColumns_, SEXP lineWeights_, SEXP keepLod_, SEXP keepLkhd_, SEXP gbLimit_, SEXP verbose_)
{
	BEGIN_RCPP
		Rcpp::NumericVector recombinationFractions;
		try
		{
			recombinationFractions = recombinationFractions_;
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input recombinationFractions must be a numeric vector");
		}
		R_xlen_t nRecombLevels = recombinationFractions.size();
		Rcpp::NumericVector::iterator halfIterator = std::find(recombinationFractions.begin(), recombinationFractions.end(), 0.5);
		if(halfIterator == recombinationFractions.end()) throw std::runtime_error("Input recombinationFractions did not contain the value 0.5");
		if(std::find(recombinationFractions.begin(), recombinationFractions.end(), 0.0) == recombinationFractions.end()) throw std::runtime_error("Input recombinationFractions did not contain the value 0");
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
		catch(...)
		{
			throw Rcpp::not_compatible("Input object must be an S4 object");
		}
		std::vector<int> markerRows;
		try
		{
			markerRows = Rcpp::as<std::vector<int> >(markerRows_);
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input markerRows must be an integer vector");
		}
		for(std::vector<int>::iterator markerRow = markerRows.begin(); markerRow != markerRows.end(); markerRow++)
		{
			(*markerRow)--;
		}

		std::vector<int> markerColumns;
		try
		{
			markerColumns = Rcpp::as<std::vector<int> >(markerColumns_);
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input markerColumns must be an integer vector");
		}
		for(std::vector<int>::iterator markerColumn = markerColumns.begin(); markerColumn != markerColumns.end(); markerColumn++)
		{
			(*markerColumn)--;
		}
		Rcpp::RObject lineWeights_noType = lineWeights_;
		if(lineWeights_noType.sexp_type() != VECSXP)
		{
			throw Rcpp::not_compatible("Input lineWeights must be a list");
		}
		Rcpp::List lineWeights = lineWeights_;
		double gbLimit;
		try
		{
			gbLimit = Rcpp::as<double>(gbLimit_);
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input gbLimit must be a single numeric value");
		}
		//The gbLimit / bytesLimit variables only apply to the generated likelihood data - They don't apply to the lookup table, which can be a huge memory bottleneck in some cases. 
		R_xlen_t bytesLimit = (R_xlen_t)(gbLimit*1000000000LL + 1LL);
	
		Rcpp::List geneticData;
		try
		{
			geneticData = object.slot("geneticData");
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input object must have a slot named \"geneticData\" which must be a list");
		}
		R_xlen_t nDesigns = geneticData.length();
		if(lineWeights.size() != nDesigns) throw std::runtime_error("Input lineWeights had the wrong number of entries");
		try
		{
			for(R_xlen_t i = 0; i < nDesigns; i++)
			{
				Rcpp::NumericVector currentDesignLineWeights = lineWeights(i);
			}
		}
		catch(...)
		{
			throw std::runtime_error("Input lineWeights must be a list of numeric vectors");
		}
		bool keepLod, keepLkhd;
		try
		{
			keepLod = Rcpp::as<bool>(keepLod_);
		}
		catch(...)
		{
			throw std::runtime_error("Input keepLod must be a boolean");
		}
		try
		{
			keepLkhd = Rcpp::as<bool>(keepLkhd_);
		}
		catch(...)
		{
			throw std::runtime_error("Input keepLkhd must be a boolean");
		}
		Rcpp::List verboseList;
		bool verbose;
		int progressStyle;
		try
		{
			verboseList = Rcpp::as<Rcpp::List>(verbose_);
			verbose = Rcpp::as<bool>(verboseList("verbose"));
			progressStyle = Rcpp::as<int>(verboseList("progressStyle"));
		}
		catch(...)
		{
			throw std::runtime_error("Input verbose must be a boolean or a list with entries verbose and progressStyle");
		}
		if (progressStyle < 1 || progressStyle > 3)
		{
			throw std::runtime_error("Input verbose$progressStyle must be 1, 2 or 3");
		}
		if(nDesigns <= 0) throw std::runtime_error("There must be at least one design");
		if(markerRows.size() == 0) throw std::runtime_error("Input markerRows must have at least one entry");
		if(markerColumns.size() == 0) throw std::runtime_error("Input markerColumns must have at least one entry");

		//Last bit of validation
		std::size_t nMarkers = 0;
		for(int i = 0; i < nDesigns; i++)
		{
			Rcpp::S4 currentGeneticData = geneticData(i);
			Rcpp::IntegerMatrix finals = currentGeneticData.slot("finals");
			if(i != 0 && (std::size_t)finals.ncol() != nMarkers)
			{
				throw std::runtime_error("Inconsistent number of markers across the designs");
			}
			nMarkers = finals.ncol();

			std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
			if((int)lineWeightsThisDesign.size() != finals.nrow())
			{
				throw std::runtime_error("An entry of input lineWeights had the wrong length");
			}
		}

		int markerRowMin = *std::min_element(markerRows.begin(), markerRows.end()), markerRowMax = *std::max_element(markerRows.begin(), markerRows.end());
		int markerColumnMin = *std::min_element(markerColumns.begin(), markerColumns.end()), markerColumnMax = *std::max_element(markerColumns.begin(), markerColumns.end());
		if(markerRowMin < 0 || (std::size_t)markerRowMax >= nMarkers) throw std::runtime_error("Invalid values for input markerRows");
		if(markerColumnMin < 0 || (std::size_t)markerColumnMax >= nMarkers) throw std::runtime_error("Invalid value for input markerColumns");

		//If the input values of markerRows and markerColumns give a region that's completely in the lower triangular region, then throw an error
		R_xlen_t nValuesToEstimate = countValuesToEstimate(markerRows, markerColumns);
		if(nValuesToEstimate == 0)
		{
			throw std::runtime_error("Input values of markerRows and markerColumns give a region that is contained in the lower triangular part of the matrix");
		}

		//Warn if we're going to allocate over 4gb
		R_xlen_t valuesToEstimateInChunk;
		if(bytesLimit < 0)
		{
			valuesToEstimateInChunk = nValuesToEstimate;
		}
		else
		{
			valuesToEstimateInChunk = std::min(nValuesToEstimate, (R_xlen_t) bytesLimit / (R_xlen_t)(nRecombLevels * sizeof(double)) + (R_xlen_t)1LL);
		}
		//Output a message detailing allocation size here, if either we're going to allocate more than 4gb, and the user didn't attempt to limit this, or the verbose option is specified. 
		if((valuesToEstimateInChunk * (std::size_t)nRecombLevels * sizeof(double) > 4000000000ULL && bytesLimit < 0) || verbose)
		{
			Rcpp::Rcout << "Allocating results matrix of " << (valuesToEstimateInChunk * (std::size_t)nRecombLevels * sizeof(double)) << " bytes = " << ((valuesToEstimateInChunk * (std::size_t)nRecombLevels * sizeof(double))/1000000000ULL) << " gb" << std::endl;
		}
		//This is not an Rcpp::NumericVector because it can quite easily overflow the size of such a vector (signed int)
		std::vector<double> result(valuesToEstimateInChunk * nRecombLevels, 0);

		//Construct vector of rfhaps_internal_args objects
		triangularIterator startPosition(markerRows, markerColumns);
		std::vector<rfhaps_internal_args> internalArgumentObjects;
		for(int i = 0; i < nDesigns; i++)
		{
			Rcpp::S4 currentGeneticData = geneticData(i);
			std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
			std::string error;
			estimateRFSpecificDesignArgs args(recombinationFractionsDouble);
			try
			{
				args.founders = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("founders"));
			}
			catch(...)
			{
				std::stringstream ss; 
				ss << "Founders slot of design " << i << " was not an integer matrix";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.finals = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("finals"));
			}
			catch(...)
			{
				std::stringstream ss; 
				ss << "Finals slot of design " << i << " was not an integer matrix";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.pedigree = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("pedigree"));
			}
			catch(...)
			{
				std::stringstream ss; 
				ss << "Pedigree slot of design " << i << " was not an S4 object";
				throw std::runtime_error(ss.str().c_str());
			}
			try
			{
				args.hetData = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("hetData"));
			}
			catch(...)
			{
				std::stringstream ss; 
				ss << "hetData slot of design " << i << " was not an S4 object";
				throw std::runtime_error(ss.str().c_str());
			}
			//This has to be copied / swapped in, because it's a local temporary at the moment
			args.lineWeights.swap(lineWeightsThisDesign);
			rfhaps_internal_args internalArgs(args.recombinationFractions, startPosition);
			internalArgs.markerRows = &markerRows;
			internalArgs.markerColumns = &markerColumns;

			bool converted = toInternalArgs(std::move(args), internalArgs, error);
			if(!converted)
			{
				std::stringstream ss;
				ss << "Error pre-processing data for dataset " << i << ": " << error;
				throw std::runtime_error(ss.str().c_str());
			}
			internalArgumentObjects.emplace_back(std::move(internalArgs));
		}
		//Estimate required memory usage
		signed long long lookupBytes = 0;
		for(int i = 0; i < nDesigns; i++)
		{
			unsigned long long currentLookupBytes = estimateLookup(internalArgumentObjects[i]);
			lookupBytes += currentLookupBytes;
		}
		//Output a message giving the allocation size, if either it's more than a gb, or the verbose option is specified
		if(lookupBytes > 1000000000 || verbose)
		{
			Rcpp::Rcout << "Total lookup table size of " << lookupBytes << " bytes = " << (lookupBytes / 1000000000ULL) << " gb" << std::endl;
		}
		Rcpp::NumericVector lod, lkhd;
		Rcpp::RawVector theta(nValuesToEstimate);
		if(keepLod) lod = Rcpp::NumericVector(nValuesToEstimate);
		if(keepLkhd) lkhd = Rcpp::NumericVector(nValuesToEstimate);
		double* resultPtr = &(result[0]);

		Rcpp::Function txtProgressBar("txtProgressBar");
		Rcpp::Function setTxtProgressBar("setTxtProgressBar");
		Rcpp::Function close("close");
		Rcpp::RObject barHandle;
		std::function<void(unsigned long long)> updateProgress = [](unsigned long long){};
		if(verbose)
		{
			barHandle = txtProgressBar(Rcpp::Named("style") = progressStyle, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
			updateProgress = [barHandle,nDesigns,nValuesToEstimate,setTxtProgressBar](unsigned long long value)
				{
					try
					{
#ifdef CUSTOM_STATIC_RCPP
						setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*value) / (double)(nDesigns*nValuesToEstimate)));
#else
						setTxtProgressBar(barHandle, (int)((double)(1000*value) / (double)(nDesigns*nValuesToEstimate)));
#endif
					}
					catch(...)
					{
					}
				};
		}
		unsigned long long counter = 0;
		for(R_xlen_t offset = 0; offset < nValuesToEstimate; offset += valuesToEstimateInChunk)
		{
			R_xlen_t valuesToEstimateInCurrentChunk = std::min(valuesToEstimateInChunk, nValuesToEstimate - offset);
			if(offset != 0) memset(resultPtr, 0, result.size() * sizeof(double));
			//Now the actual computation)
			for(int i = 0; i < nDesigns; i++)
			{
				internalArgumentObjects[i].result = resultPtr;
				internalArgumentObjects[i].valuesToEstimateInChunk = valuesToEstimateInCurrentChunk;
				internalArgumentObjects[i].startPosition = startPosition;
				internalArgumentObjects[i].updateProgress = updateProgress;
				std::string error;
				bool successful = estimateRFSpecificDesign(internalArgumentObjects[i], counter);
				if(!successful) throw std::runtime_error("Internal error");
			}
			//now for some post-processing to get out the MLE, lod (maybe) and lkhd (maybe)
			R_xlen_t endValue = std::min(nValuesToEstimate, offset + valuesToEstimateInCurrentChunk);
			for(R_xlen_t counter = offset; counter < endValue; counter++)
			{
				double* start = resultPtr + (counter - offset) * nRecombLevels;
				double* end = resultPtr + (counter - offset + 1) * nRecombLevels;
				double* maxPtr = std::max_element(start, end), *minPtr = std::min_element(start, end);
				double max = *maxPtr, min = *minPtr;
				int currentTheta;
				double currentLod;
				//This is the case where no data was available, across any of the experiments. This is precise, no numerical error involved
				if(max == 0 && min == 0)
				{
					max = currentLod = std::numeric_limits<double>::quiet_NaN();
					currentTheta = 0xff;
				}
				else
				{
					currentTheta = (int)(maxPtr - start);
					currentLod = max - resultPtr[(counter - offset)* (std::ptrdiff_t)nRecombLevels + halfIndex];
				}
				theta(counter) = currentTheta;
				if(keepLkhd) lkhd(counter) = max;
				if(keepLod) lod(counter) = currentLod;
			}
			while(valuesToEstimateInCurrentChunk > 0)
			{
				startPosition.next();
				valuesToEstimateInCurrentChunk--;
			}
		}
		if(verbose)
		{
			close(barHandle);
		}
		Rcpp::RObject lodRet, lkhdRet;
		
		if(keepLod) lodRet = lod;
		else lodRet = R_NilValue;

		if(keepLkhd) lkhdRet = lkhd;
		else lkhdRet = R_NilValue;

		return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lod") = lodRet, Rcpp::Named("lkhd") = lkhdRet, Rcpp::Named("r") = recombinationFractions);
	END_RCPP
}
