#include "estimateRF.h"
#include <math.h>
#include <limits>
#include "estimateRFSingleDesign.h"
#include "estimateRFSingleDesignInternal.h"
#include <stdexcept>
#include "matrixChunks.h"
SEXP estimateRFSingleDesign(SEXP object_, SEXP recombinationFractions_, SEXP markerRows_, SEXP markerColumns_, SEXP lineWeights_, SEXP keepLod_, SEXP keepLkhd_, SEXP verbose_)
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
		Rcpp::NumericVector::iterator halfIterator = std::find(recombinationFractions.begin(), recombinationFractions.end(), 0.5);
		if(halfIterator == recombinationFractions.end()) throw std::runtime_error("Input recombinationFractions did not contain the value 0.5");
		if(std::find(recombinationFractions.begin(), recombinationFractions.end(), 0.0) == recombinationFractions.end()) throw std::runtime_error("Input recombinationFractions did not contain the value 0");

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
		//Change from 1-base indexing to 0-base indexing. 
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
		//Change from 1-base indexing to 0-base indexing. 
		for(std::vector<int>::iterator markerColumn = markerColumns.begin(); markerColumn != markerColumns.end(); markerColumn++)
		{
			(*markerColumn)--;
		}

		Rcpp::List geneticData;
		try
		{
			geneticData = object.slot("geneticData");
		}
		catch(...)
		{
			throw Rcpp::not_compatible("Input object must have a slot named \"geneticData\" which must be a list");
		}
		if(geneticData.size() != 1)
		{
			throw std::runtime_error("Only a single genetic design can be input");
		}

		Rcpp::NumericVector lineWeights;
		try
		{
			lineWeights = lineWeights_;
		}
		catch(...)
		{
			throw std::runtime_error("Input lineWeights must be a numeric vector");
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
		if(markerRows.size() == 0) throw std::runtime_error("Input markerRows must have at least one entry");
		if(markerColumns.size() == 0) throw std::runtime_error("Input markerColumns must have at least one entry");

		Rcpp::S4 currentGeneticData = geneticData(0);
		Rcpp::IntegerMatrix finals = currentGeneticData.slot("finals");
		std::size_t nMarkers = finals.ncol();

		std::vector<double> lineWeightsDouble = Rcpp::as<std::vector<double> >(lineWeights);
		if((int)lineWeightsDouble.size() != finals.nrow())
		{
			throw std::runtime_error("An entry of input lineWeights had the wrong length");
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

		std::string error;
		estimateRFSingleDesignArgs args(recombinationFractionsDouble);
		try
		{
			args.founders = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("founders"));
		}
		catch(...)
		{
			std::stringstream ss; 
			ss << "Founders slot was not an integer matrix";
			throw std::runtime_error(ss.str().c_str());
		}
		try
		{
			args.finals = Rcpp::as<Rcpp::IntegerMatrix>(currentGeneticData.slot("finals"));
		}
		catch(...)
		{
			std::stringstream ss; 
			ss << "Finals slot was not an integer matrix";
			throw std::runtime_error(ss.str().c_str());
		}
		try
		{
			args.pedigree = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("pedigree"));
		}
		catch(...)
		{
			std::stringstream ss; 
			ss << "Pedigree slot was not an S4 object";
			throw std::runtime_error(ss.str().c_str());
		}
		try
		{
			args.hetData = Rcpp::as<Rcpp::S4>(currentGeneticData.slot("hetData"));
		}
		catch(...)
		{
			std::stringstream ss; 
			ss << "hetData slot was not an S4 object";
			throw std::runtime_error(ss.str().c_str());
		}
		//This has to be copied / swapped in, because it's a local temporary at the moment
		args.lineWeights.swap(lineWeightsDouble);
		estimateRFSingleDesignInternalArgs internalArgs(args.recombinationFractions);
		internalArgs.markerRows = &markerRows;
		internalArgs.markerColumns = &markerColumns;

		bool converted = toInternalArgs(std::move(args), internalArgs, error);
		if(!converted)
		{
			std::stringstream ss;
			ss << "Error pre-processing dataset: " << error;
			throw std::runtime_error(ss.str().c_str());
		}

		Rcpp::NumericVector lod, lkhd;
		Rcpp::RawVector theta(nValuesToEstimate);
		if(keepLod) 
		{
			lod = Rcpp::NumericVector(nValuesToEstimate);
			internalArgs.lod = &(lod[0]);
		}
		if(keepLkhd)
		{
			lkhd = Rcpp::NumericVector(nValuesToEstimate);
			internalArgs.lkhd = &(lkhd[0]);
		}
		internalArgs.theta = &(theta[0]);

		Rcpp::Environment progress = Rcpp::Environment::namespace_env("progress");
		Rcpp::Environment progressBar = progress.get("progress_bar");
		Rcpp::Function newProgressBar = progressBar["new"];
		Rcpp::Environment barHandle;
		std::function<void(unsigned long long)> updateProgress = [](unsigned long long){};
		if(verbose)
		{
			barHandle = newProgressBar(Rcpp::Named("format") = "Estimating [:bar] :percent eta :eta", Rcpp::Named("total") = (double)nValuesToEstimate, Rcpp::Named("clear") = true);
			Rcpp::Function tick = barHandle["tick"];
			unsigned long long previousValue = 0;
			updateProgress = [barHandle, tick, &previousValue](unsigned long long value)
				{
					try
					{
						int increment = (int)(value - previousValue);
#ifdef CUSTOM_STATIC_RCPP
						tick.topLevelExec(increment);
#else
						tick(increment);
#endif
						previousValue = value;
					}
					catch(...)
					{
					}
				};
		}
		internalArgs.updateProgress = updateProgress;
		estimateRFSingleDesignInternal(internalArgs);
		Rcpp::RObject lodRet, lkhdRet;
		
		if(keepLod) lodRet = lod;
		else lodRet = R_NilValue;

		if(keepLkhd) lkhdRet = lkhd;
		else lkhdRet = R_NilValue;

		return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lod") = lodRet, Rcpp::Named("lkhd") = lkhdRet, Rcpp::Named("r") = recombinationFractions);
	END_RCPP
}
