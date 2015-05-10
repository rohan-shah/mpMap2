#include "estimateRf.h"
#include <math.h>
#include <limits>
#include "estimateRfSpecificDesign.h"
#include <stdexcept>
SEXP estimateRf(SEXP object_, SEXP recombinationFractions_, SEXP marker1Range_, SEXP marker2Range_, SEXP lineWeights_, SEXP keepLod_, SEXP keepLkhd_)
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
		int halfIndex = std::distance(recombinationFractions.begin(), halfIterator);

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
		Rcpp::IntegerVector marker1Range;
		try
		{
			marker1Range = marker1Range_;
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input marker1Range must be an integer vector");
		}
		Rcpp::IntegerVector marker2Range;
		try
		{
			marker2Range = marker2Range_;
		}
		catch(Rcpp::not_compatible&)
		{
			throw Rcpp::not_compatible("Input marker2Range must be an integer vector");
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
		if(marker1Range.size() != 2) throw std::runtime_error("Input marker1Range must have two entries");
		if(marker2Range.size() != 2) throw std::runtime_error("Input marker2Range must have two entries");

		int marker1Start = marker1Range(0), marker1End = marker1Range(1);
		int marker2Start = marker2Range(0), marker2End = marker2Range(1);
		if(marker1Start > marker1End || marker1Start < 0 || marker1End < 0) throw std::runtime_error("Invalid value for input marker1Range");
		if(marker2Start > marker2End || marker2Start < 0 || marker2End < 0) throw std::runtime_error("Invalid value for input marker2Range");
		long marker1RangeSize = marker1End - marker1Start, marker2RangeSize = marker2End - marker2Start;

		//Indexing has form result[markerCounter1 *nRecombLevels*nMarkers + markerCounter2 * nRecombLevels + recombCounter]
		//This is not an Rcpp::NumericVector because it can quite easily overflow the size of such a vector (signed int)
		std::vector<double> result(marker1RangeSize * marker2RangeSize * nRecombLevels, 0);

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
		Rcpp::CharacterVector allMarkerNames = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::RObject>(Rcpp::as<Rcpp::S4>(geneticData(0)).slot("finals")).attr("dimnames"))[1];
		//Now the actual computation
		for(int i = 0; i < nDesigns; i++)
		{
			Rcpp::S4 currentGeneticData = geneticData(i);
			std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
			std::string error;
			estimateRfSpecificDesignArgs args(lineWeightsThisDesign, recombinationFractionsDouble);
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
			args.marker1Start = marker1Start;
			args.marker2Start = marker2Start;
			args.marker1End = marker1End;
			args.marker2End = marker2End;
			args.result = &(result[0]);
			bool successful = estimateRfSpecificDesign(args);
			if(!successful) throw std::runtime_error(args.error);
		}
		//now for some post-processing to get out the MLE, lod (maybe) and lkhd (maybe)
		Rcpp::NumericMatrix theta(marker1RangeSize, marker2RangeSize), lod, lkhd;
		if(keepLod) lod = Rcpp::NumericMatrix(marker1RangeSize, marker2RangeSize);
		if(keepLkhd) lkhd = Rcpp::NumericMatrix(marker1RangeSize, marker2RangeSize);
	
		//row and column names of output
		Rcpp::CharacterVector outputRowNames(marker1RangeSize), outputColumnNames(marker2RangeSize);
		std::copy(allMarkerNames.begin() + marker1Start, allMarkerNames.begin() + marker1Start + marker1RangeSize, outputRowNames.begin());
		std::copy(allMarkerNames.begin() + marker2Start, allMarkerNames.begin() + marker2Start + marker2RangeSize, outputColumnNames.begin());

		Rcpp::List outputDimNames(Rcpp::List::create(outputRowNames, outputColumnNames));
		theta.attr("dimnames") = outputDimNames;
		if(keepLod) lod.attr("dimnames") = outputDimNames;
		if(keepLkhd) lkhd.attr("dimnames") = outputDimNames;

		double* resultPtr = &(result[0]);
		int maxStart = std::max(marker1Start, marker2Start);
		int minEnd = std::min(marker1End, marker2End);
		//first get maximum and minimum likelihoods
		#ifdef USE_OPENMP
			#pragma omp parallel for
		#endif
		for(int markerCounter1 = marker1Start; markerCounter1 < marker1End; markerCounter1++)
		{
			for(int markerCounter2 = marker2Start; markerCounter2 < marker2End; markerCounter2++)
			{
				int destIndex1 = markerCounter1 - marker1Start, destIndex2 = markerCounter2 - marker2Start;
				int sourceIndex1 = destIndex1, sourceIndex2 = destIndex2;
				if(markerCounter2 >= maxStart && markerCounter2 < minEnd && markerCounter1 >= maxStart && markerCounter1 < minEnd && markerCounter2 < markerCounter1)
				{
					sourceIndex1 = markerCounter2 - maxStart + (maxStart - marker1Start);
					sourceIndex2 = markerCounter1 - maxStart + (maxStart - marker2Start);
				}

				double* start = resultPtr + (long)sourceIndex1 *(long)marker2RangeSize*(long)nRecombLevels + (long)sourceIndex2 * (long)nRecombLevels;
				double* end = resultPtr + (long)sourceIndex1 *(long)marker2RangeSize*(long)nRecombLevels + (long)sourceIndex2 * (long)nRecombLevels + (long)nRecombLevels;
				double* maxPtr = std::max_element(start, end), *minPtr = std::min_element(start, end);
				double max = *maxPtr, min = *minPtr;
				double currentTheta, currentLod;
				//This is the case where no data was available, across any of the experiments. This is precise, no numerical error involved
				if(max == 0 && min == 0)
				{
					currentLod = currentTheta = std::numeric_limits<double>::quiet_NaN();
				}
				else
				{
					currentTheta = recombinationFractions(maxPtr - start);
					currentLod = max - resultPtr[(long)sourceIndex1*(long)nRecombLevels*(long)marker2RangeSize + (long)sourceIndex2 * (long)nRecombLevels + halfIndex];
				}
				theta(destIndex1, destIndex2) = currentTheta;
				if(keepLkhd) lkhd(destIndex1, destIndex2) = max;
				if(keepLod) lod(destIndex1, destIndex2) = currentLod;
			}
		}
		Rcpp::RObject lodRet, lkhdRet;
		
		if(keepLod) lodRet = lod;
		else lodRet = R_NilValue;

		if(keepLkhd) lkhdRet = lkhd;
		else lkhdRet = R_NilValue;

		return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lod") = lodRet, Rcpp::Named("lkhd") = lkhdRet, Rcpp::Named("r") = recombinationFractions);
	END_RCPP
}
