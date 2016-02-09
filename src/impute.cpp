#include "impute.h"
#include <vector>
#include <map>
#include <math.h>
#include <limits>
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif
template<bool hasLOD, bool hasLKHD> bool imputeInternal(unsigned char* theta, std::vector<double>& levels, double* lod, double* lkhd, std::vector<int>& markersThisGroup, std::string& error, std::function<void(unsigned long, unsigned long)> statusFunction)
{
	unsigned long done = 0;
	unsigned long total = markersThisGroup.size();
	unsigned long doneThreadZero = 0;
	bool hasError = false;
	//This is a marker row
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for(std::vector<int>::iterator marker1 = markersThisGroup.begin(); marker1 != markersThisGroup.end(); marker1++)
	{
		bool missing = false;
		for(std::vector<int>::iterator marker2 = markersThisGroup.begin(); marker2 != markersThisGroup.end(); marker2++)
		{
			int copiedMarker1 = *marker1;
			int copiedMarker2 = *marker2;
			if(copiedMarker1 > copiedMarker2) std::swap(copiedMarker1, copiedMarker2);
			//record whether there's a missing value for marker1 anywhere.
			float val = theta[(copiedMarker2 * (copiedMarker2 + 1))/2 + copiedMarker1];
			if(val == 0xff) 
			{
				missing = true;
				break;
			}
		}
		if(missing)
		{
			//There's a missing value for *marker1. So find the closest matching column and overwrite any missing values from that column.
			//Map with the difference as the key and the marker index as the value. These are sorted internally so we can go from smallest difference to lowest difference
			std::map<float, int> averageDifferences;
			//marker 2 is the candidate other marker (another row)
			for(std::vector<int>::iterator marker2 = markersThisGroup.begin(); marker2 != markersThisGroup.end(); marker2++)
			{
				if(marker2 == marker1) continue;
				
				float totalDifference = 0;
				int usableLocations = 0;
				//This iterates over the other markers in the column/row
				for(std::vector<int>::iterator marker3 = markersThisGroup.begin(); marker3 != markersThisGroup.end(); marker3++)
				{
					int pair1Row = *marker1, pair1Column = *marker3;
					int pair2Row = *marker2, pair2Column = *marker3;
					if(pair1Row > pair1Column) std::swap(pair1Row, pair1Column);
					if(pair2Row > pair2Column) std::swap(pair2Row, pair2Column);
					unsigned char value1 = theta[(pair1Column * (pair1Column + 1))/2 + pair1Row];
					unsigned char value2 = theta[(pair2Column * (pair2Column + 1))/2 + pair2Row];
					if(value1 != 0xff && value2 != 0xff)
					{
						float val = (float)fabs(levels[value1] - levels[value2]);
						totalDifference += val;
						usableLocations++;
					}
				}
				if(usableLocations != 0)
				{
					averageDifferences.insert(std::make_pair(totalDifference/usableLocations, *marker2));
				}
				else averageDifferences.insert(std::make_pair(std::numeric_limits<float>::quiet_NaN(), *marker2));
			}
			for(std::vector<int>::iterator marker2 = markersThisGroup.begin(); marker2 != markersThisGroup.end(); marker2++)
			{
				int pair1Row = *marker1;
				int pair1Column = *marker2;
				if(pair1Row > pair1Column) std::swap(pair1Row, pair1Column);
				unsigned char& toReplace = theta[(pair1Column*(pair1Column + 1))/2 + pair1Row];
				//if this is missing, we need to impute it
				if(toReplace == 0xff) 
				{
					//go through the other markers from most similar to least similar, looking for something which has a value here. So marker3 is the 
					bool replacementFound = false;
					for(std::map<float, int>::iterator marker3 = averageDifferences.begin(); marker3 != averageDifferences.end(); marker3++)
					{
						int pair2Row = marker3->second;
						int pair2Column = *marker2;
						if(pair2Row > pair2Column) std::swap(pair2Row, pair2Column);
						unsigned char& newValue = theta[(pair2Column * (pair2Column + 1))/2 + pair2Row];
						if(newValue != 0xff)
						{
							toReplace = newValue;
							if(hasLOD) lod[(pair1Column*(pair1Column + 1))/2 + pair1Row] = lod[(pair2Column * (pair2Column + 1))/2 + pair2Row];
							if(hasLKHD) lkhd[(pair1Column*(pair1Column + 1))/2 + pair1Row] = lkhd[(pair2Column * (pair2Column + 1))/2 + pair2Row];
							//We only need to copy the value from the best other similar marker, so we can break here
							replacementFound = true;
							break;
						}
					}
					if(!replacementFound)
					{
						std::stringstream ss;
						ss << "Unable to impute a value for marker " << (*marker1+1) << " and marker " << (*marker2+1);
						error = ss.str();
						hasError = true;
#ifndef USE_OPENMP
						//We can't return early if we're using openmp
						return false;
#endif
					}
				}
			}
		}
#ifdef USE_OPENMP
		#pragma omp critical
#endif
		{
			done++;
		}
#ifdef USE_OPENMP
		if(omp_get_thread_num() == 0)
#endif
		{
			doneThreadZero++;
			statusFunction(done, total);
		}
	}
	return !hasError;
}
bool impute(unsigned char* theta, std::vector<double>& thetaLevels, double* lod, double* lkhd, std::vector<int>& markers, std::string& error, std::function<void(unsigned long, unsigned long)> statusFunction)
{
	if(lod != NULL && lkhd != NULL)
	{
		return imputeInternal<true, true>(theta, thetaLevels, lod, lkhd, markers, error, statusFunction);
	}
	else if(lod != NULL && lkhd == NULL)
	{
		return imputeInternal<true, false>(theta, thetaLevels, lod, lkhd, markers, error, statusFunction);
	}
	else if(lod == NULL && lkhd != NULL)
	{
		return imputeInternal<false, true>(theta, thetaLevels, lod, lkhd, markers, error, statusFunction);
	}
	else
	{
		return imputeInternal<false, false>(theta, thetaLevels, lod, lkhd, markers, error, statusFunction);
	}
}
SEXP imputeWholeObject(SEXP mpcrossLG_sexp, SEXP verbose_sexp)
{
BEGIN_RCPP
	Rcpp::S4 mpcrossLG;
	try
	{
		mpcrossLG = Rcpp::as<Rcpp::S4>(mpcrossLG_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input mpcrossLG must be an S4 object");
	}

	Rcpp::S4 rf;
	try
	{
		rf = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("rf"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf must be an S4 object");
	}

	Rcpp::S4 theta;
	try
	{
		theta = Rcpp::as<Rcpp::S4>(rf.slot("theta"));
	}
	catch(...)
	{
		 throw std::runtime_error("Slot mpcrossLG@rf@theta must be an S4 object");
	}

	Rcpp::RawVector thetaData;
	try
	{
		thetaData = Rcpp::as<Rcpp::RawVector>(theta.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@data must be a raw vector");
	}
	
	Rcpp::RawVector copiedThetaData(thetaData.size());
	memcpy(&(copiedThetaData[0]), &(thetaData[0]), sizeof(Rbyte)*thetaData.size());

	std::vector<double> levels;
	try
	{
		levels = Rcpp::as<std::vector<double> >(theta.slot("levels"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@levels must be an integer vector");
	}

	Rcpp::NumericVector copiedLod, copiedLkhd;
	if(!Rcpp::as<Rcpp::RObject>(rf.slot("lkhd")).isNULL())
	{
		Rcpp::S4 lkhdS4;
		try
		{
			lkhdS4 = Rcpp::as<Rcpp::S4>(rf.slot("lkhd"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lkhd must be an S4 object or NULL");
		}
		
		Rcpp::NumericVector lkhdS4Data;
		try
		{
			lkhdS4Data = Rcpp::as<Rcpp::NumericVector>(lkhdS4.slot("x"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lkhd@x must be a numeric vector");
		}
		copiedLkhd = Rcpp::NumericVector(lkhdS4Data.size());
		memcpy(&(copiedLkhd[0]), &(lkhdS4Data[0]), sizeof(double)*lkhdS4Data.size());
	}

	if(!Rcpp::as<Rcpp::RObject>(rf.slot("lod")).isNULL())
	{
		Rcpp::S4 lodS4;
		try
		{
			lodS4 = Rcpp::as<Rcpp::S4>(rf.slot("lod"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lod must be an S4 object or NULL");
		}
		
		Rcpp::NumericVector lodS4Data;
		try
		{
			lodS4Data = Rcpp::as<Rcpp::NumericVector>(lodS4.slot("x"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lod@x must be a numeric vector");
		}
		copiedLod = Rcpp::NumericVector(lodS4Data.size());
		memcpy(&(copiedLod[0]), &(lodS4Data[0]), sizeof(double)*lodS4Data.size());
	}

	Rcpp::S4 lg;
	try
	{
		lg = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("lg"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcross@lg must be an S4 object");
	}

	Rcpp::IntegerVector groups;
	std::vector<int> allGroups;
	try
	{
		groups = Rcpp::as<Rcpp::IntegerVector>(lg.slot("groups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcross@lg@groups must be an integer vector");
	}

	try
	{
		allGroups = Rcpp::as<std::vector<int> >(lg.slot("allGroups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcross@lg@allGroups must be an integer vector");
	}

	bool verbose;
	try
	{
		verbose = Rcpp::as<bool>(verbose_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input verbose must be a boolean");
	}

	std::vector<int> markersCurrentGroup;
	double* lodPtr = NULL, *lkhdPtr = NULL;
	if(copiedLod.size() != 0)
	{
		lodPtr = &(copiedLod[0]);
	}
	if(copiedLkhd.size() != 0)
	{
		lkhdPtr = &(copiedLkhd[0]);
	}

	Rcpp::Function txtProgressBar("txtProgressBar"), setTxtProgressBar("setTxtProgressBar"), close("close");
	Rcpp::RObject barHandle;
	std::function<void(unsigned long, unsigned long)> progressFunction = [](unsigned long, unsigned long){};
	for(std::vector<int>::iterator group = allGroups.begin(); group != allGroups.end(); group++)
	{
		markersCurrentGroup.clear();
		for(R_xlen_t markerCounter = 0; markerCounter < groups.size(); markerCounter++)
		{
			if(groups[markerCounter] == *group)
			{
				markersCurrentGroup.push_back((int)markerCounter);
			}
		}
		if(verbose)
		{
			Rcpp::Rcout << "Starting imputation for group " << *group << std::endl;
			barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
			progressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
			{
				setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
			};
		}

		std::string error;
		bool ok = impute(&(copiedThetaData[0]), levels, lodPtr, lkhdPtr, markersCurrentGroup, error, progressFunction);
		if(!ok)
		{
			std::stringstream ss;
			ss << "Error performing imputation for group " << *group << ": " << error;
			throw std::runtime_error(ss.str().c_str());
		}
		if(verbose)
		{
			close(barHandle);
		}
	}
	return Rcpp::List::create(Rcpp::Named("theta") = copiedThetaData, Rcpp::Named("lod") = copiedLod, Rcpp::Named("lkhd") = copiedLkhd);
END_RCPP
}
SEXP imputeGroup(SEXP mpcrossLG_sexp, SEXP verbose_sexp, SEXP group_sexp)
{
BEGIN_RCPP
	Rcpp::S4 mpcrossLG;
	try
	{
		mpcrossLG = Rcpp::as<Rcpp::S4>(mpcrossLG_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input mpcrossLG must be an S4 object");
	}

	Rcpp::S4 rf;
	try
	{
		rf = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("rf"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf must be an S4 object");
	}

	Rcpp::S4 theta;
	try
	{
		theta = Rcpp::as<Rcpp::S4>(rf.slot("theta"));
	}
	catch(...)
	{
		 throw std::runtime_error("Slot mpcrossLG@rf@theta must be an S4 object");
	}

	Rcpp::RawVector thetaData;
	try
	{
		thetaData = Rcpp::as<Rcpp::RawVector>(theta.slot("data"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@data must be a raw vector");
	}
	
	std::vector<double> levels;
	try
	{
		levels = Rcpp::as<std::vector<double> >(theta.slot("levels"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@levels must be an integer vector");
	}

	Rcpp::NumericVector copiedLod, copiedLkhd;
	double *copiedLodPtr = NULL, *copiedLkhdPtr = NULL;
	Rcpp::NumericVector lodS4Data, lkhdS4Data;
	if(!Rcpp::as<Rcpp::RObject>(rf.slot("lkhd")).isNULL())
	{
		Rcpp::S4 lkhdS4;
		try
		{
			lkhdS4 = Rcpp::as<Rcpp::S4>(rf.slot("lkhd"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lkhd must be an S4 object or NULL");
		}
		
		try
		{
			lkhdS4Data = Rcpp::as<Rcpp::NumericVector>(lkhdS4.slot("x"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lkhd@x must be a numeric vector");
		}
		copiedLkhd = Rcpp::NumericVector(lkhdS4Data.size());
		copiedLkhdPtr = &(copiedLkhd[0]);
	}

	if(!Rcpp::as<Rcpp::RObject>(rf.slot("lod")).isNULL())
	{
		Rcpp::S4 lodS4;
		try
		{
			lodS4 = Rcpp::as<Rcpp::S4>(rf.slot("lod"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lod must be an S4 object or NULL");
		}
		
		try
		{
			lodS4Data = Rcpp::as<Rcpp::NumericVector>(lodS4.slot("x"));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@lod@x must be a numeric vector");
		}
		copiedLod = Rcpp::NumericVector(lodS4Data.size());
		copiedLodPtr = &(copiedLod[0]);
	}

	Rcpp::S4 lg;
	try
	{
		lg = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("lg"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcross@lg must be an S4 object");
	}

	Rcpp::IntegerVector groups;
	try
	{
		groups = Rcpp::as<Rcpp::IntegerVector>(lg.slot("groups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcross@lg@groups must be an integer vector");
	}

	int group;
	try
	{
		group = Rcpp::as<int>(group_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input group must be an integer");
	}

	std::vector<int> markersCurrentGroup;
	for(R_xlen_t markerCounter = 0; markerCounter < groups.size(); markerCounter++)
	{
		if(groups[markerCounter] == group)
		{
			markersCurrentGroup.push_back((int)markerCounter);
		}
	}
	if(markersCurrentGroup.size() == 0)
	{
		throw std::runtime_error("No markers belonged to the specified group");
	}

	bool verbose;
	try
	{
		verbose = Rcpp::as<bool>(verbose_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input verbose must be a boolean");
	}

	Rcpp::RawVector copiedTheta((markersCurrentGroup.size()*(markersCurrentGroup.size() + 1))/2);
	//column
	for(int marker1Counter = 0; marker1Counter < markersCurrentGroup.size(); marker1Counter++)
	{
		//row
		for(int marker2Counter = 0; marker2Counter <= marker1Counter; marker2Counter++)
		{
			copiedTheta[(marker1Counter*(marker1Counter+1))/2 + marker2Counter] = thetaData[(markersCurrentGroup[marker1Counter] *(markersCurrentGroup[marker1Counter] + 1))/2 + markersCurrentGroup[marker2Counter]];
			if(copiedLodPtr) copiedLodPtr[(marker1Counter*(marker1Counter+1))/2 + marker2Counter] = lodS4Data[(markersCurrentGroup[marker1Counter] *(markersCurrentGroup[marker1Counter] + 1))/2 + markersCurrentGroup[marker2Counter]];
			if(copiedLkhdPtr) copiedLkhdPtr[(marker1Counter*(marker1Counter+1))/2 + marker2Counter] = lkhdS4Data[(markersCurrentGroup[marker1Counter] *(markersCurrentGroup[marker1Counter] + 1))/2 + markersCurrentGroup[marker2Counter]];
		}
	}

	std::function<void(unsigned long, unsigned long)> progressFunction = [](unsigned long, unsigned long){};
	Rcpp::Function txtProgressBar("txtProgressBar"), setTxtProgressBar("setTxtProgressBar"), close("close");
	Rcpp::RObject barHandle;
	if(verbose)
	{
		Rcpp::Rcout << "Starting imputation for group " << group << std::endl;
		barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
		progressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
		{
			setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
		};
	}
	std::string error;
	bool ok = impute(&(copiedTheta[0]), levels, copiedLodPtr, copiedLkhdPtr, markersCurrentGroup, error, progressFunction);
	if(!ok)
	{
		std::stringstream ss;
		ss << "Error performing imputation for group " << group << ": " << error;
		throw std::runtime_error(ss.str().c_str());
	}
	if(verbose)
	{
		close(barHandle);
	}
	return Rcpp::List::create(Rcpp::Named("theta") = copiedTheta, Rcpp::Named("lod") = copiedLod, Rcpp::Named("lkhd") = copiedLkhd);
END_RCPP
}
