#include "order.h"
#include "impute.h"
#include "arsaRaw.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
SEXP order(SEXP mpcrossLG_sexp, SEXP groupsToOrder_sexp, SEXP cool_, SEXP temperatureMin_, SEXP nReps_, SEXP maxMove_sexp, SEXP effortMultiplier_sexp, SEXP randomStart_sexp, SEXP verbose_)
{
BEGIN_RCPP
	Rcpp::S4 mpcrossLG;
	try
	{
		mpcrossLG = mpcrossLG_sexp;
	}
	catch(...)
	{
		throw std::runtime_error("Input mpcrossLG must be an S4 object");
	}
	Rcpp::Function validObject("validObject");
	try
	{
		validObject(mpcrossLG);
	}
	catch(...)
	{
		throw std::runtime_error("Input mpcrossLG object was invalid. Please call validObject for more details");
	}

	Rcpp::S4 lg;
	try
	{
		lg = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("lg"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg must be an S4 object");
	}

	bool randomStart;
	try
	{
		randomStart = Rcpp::as<bool>(randomStart_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input randomStart must be a logical");
	}

	int maxMove;
	try
	{
		maxMove = Rcpp::as<int>(maxMove_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input maxMove must be an integer");
	}
	if(maxMove < 0)
	{
		throw std::runtime_error("Input maxMove must be non-negative");
	}

	double effortMultiplier;
	try
	{
		effortMultiplier = Rcpp::as<double>(effortMultiplier_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input effortMultiplier must be numeric");
	}
	if(effortMultiplier <= 0)
	{
		throw std::runtime_error("Input effortMultiplier must be positive");
	}

	std::vector<int> groups;
	try
	{
		groups = Rcpp::as<std::vector<int> >(lg.slot("groups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg@groups must be an integer vector");
	}

	std::vector<int> allGroups;
	try
	{
		allGroups = Rcpp::as<std::vector<int> >(lg.slot("allGroups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg@allGroups must be an integer vector");
	}

	std::vector<int> groupsToOrder;
	try
	{
		groupsToOrder = Rcpp::as<std::vector<int> >(groupsToOrder_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input groupsToOrder must be an integerVector");
	}

	double cool;
	try
	{
		cool = Rcpp::as<double>(cool_);
	}
	catch(...)
	{
		throw std::runtime_error("Input cool must be a number");
	}

	double temperatureMin;
	try
	{
		temperatureMin = Rcpp::as<double>(temperatureMin_);
	}
	catch(...)
	{
		throw std::runtime_error("Input temperatureMin must be a number");
	}

	long nReps;
	try
	{
		nReps = Rcpp::as<int>(nReps_);
	}
	catch(...)
	{
		throw std::runtime_error("Input nReps must be an integer");
	}

	bool verbose;
	try
	{
		verbose = Rcpp::as<bool>(verbose_);
	}
	catch(...)
	{
		throw std::runtime_error("Input verbose must be a boolean");
	}
	std::vector<double> levels;


	Rcpp::RObject imputedTheta_robject;
	try
	{
		imputedTheta_robject = lg.slot("imputedTheta");
	}
	catch(...)
	{
		throw std::runtime_error("Internal error accessing slot mpcrossLG@lg@imputedTheta");
	}
	bool hasImputedTheta = !imputedTheta_robject.isNULL();
	Rcpp::RawVector thetaRawData;
	Rcpp::List imputedTheta;
	if(hasImputedTheta)
	{
		imputedTheta = Rcpp::as<Rcpp::List>(imputedTheta_robject);
		if(verbose)
		{
			Rcpp::Rcout << "Using previously imputed recombination fraction matrices" << std::endl;
		}
		if((int)imputedTheta.size() != (int)allGroups.size())
		{
			throw std::runtime_error("Slot mpcrossLG@lg@imputedTheta had the wrong length");
		}
		if(imputedTheta.size() > 0)
		{
			levels = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::S4>(imputedTheta(0)).slot("levels"));
		}
	}
	else
	{
		Rcpp::S4 rf;
		try
		{
			rf = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("rf"));
		}
		catch(...)
		{
			throw std::runtime_error("If slot mpcrossLG@lg@imputedTheta is missing, mpcrossLG@rf cannot be missing and must be an S4 object");
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

		try
		{
			thetaRawData = theta.slot("data");
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@theta@data must be a raw vector");
		}

		try
		{
			levels = Rcpp::as<std::vector<double> >(Rcpp::as<Rcpp::NumericVector>(theta.slot("levels")));
		}
		catch(...)
		{
			throw std::runtime_error("Slot mpcrossLG@rf@theta@levels must be a numeric vector");
		}
	}

	//Check that groupsToOrder contains unique values, which are contained in allGroups
	std::sort(groupsToOrder.begin(), groupsToOrder.end());
	std::sort(allGroups.begin(), allGroups.end());
	if(std::unique(groupsToOrder.begin(), groupsToOrder.end()) != groupsToOrder.end())
	{
		throw std::runtime_error("Input groupsToOrder must contain unique values");
	}
	for(std::vector<int>::iterator groupToOrder = groupsToOrder.begin(); groupToOrder != groupsToOrder.end(); groupToOrder++)
	{
		std::vector<int>::iterator bound = std::lower_bound(allGroups.begin(), allGroups.end(), *groupToOrder);
		if(bound == allGroups.end() || *bound != *groupToOrder)
		{
			throw std::runtime_error("Input groupsToOrder contained a group that was not present in slot mpcrossLG@lg@allGroups");
		}
	}


	R_xlen_t nMarkers = groups.size();

	std::vector<int> permutation, currentGroupPermutation;
	permutation.reserve(nMarkers);
	std::vector<int> markersThisGroup;
	//This holds a copy of the raw data, for the purposes of doing imputation. We don't want to touch the original, obviously. 
	std::vector<unsigned char> imputedRaw;
	unsigned char* imputedRawPtr;
	//Stuff for the verbose output case
	Rcpp::RObject barHandle;
	Rcpp::Function txtProgressBar("txtProgressBar"), setTxtProgressBar("setTxtProgressBar"), close("close");
	for(std::vector<int>::iterator currentGroup = allGroups.begin(); currentGroup != allGroups.end(); currentGroup++)
	{
		int groupCount = (int)std::distance(allGroups.begin(), currentGroup);
		markersThisGroup.clear();
		//linkage groups are not required to be contiguous in this case, so we have to scan through to find the number of markers in this group
		for(std::vector<int>::iterator currentMarkerGroup = groups.begin(); currentMarkerGroup != groups.end(); currentMarkerGroup++)
		{
			if(*currentMarkerGroup == *currentGroup)
			{
				markersThisGroup.push_back((int)std::distance(groups.begin(), currentMarkerGroup));
			}
		}
		std::size_t nMarkersCurrentGroup = (int)markersThisGroup.size();
		if(nMarkersCurrentGroup == 0) continue;
		
		std::vector<int> contiguousIndices(nMarkersCurrentGroup);
		for(std::size_t i = 0; i < nMarkersCurrentGroup; i++) contiguousIndices[i] = (int)i;

		if(!hasImputedTheta)
		{
			//Do imputation
			//So first make a copy of the raw data subset
			imputedRaw.resize((nMarkersCurrentGroup * (nMarkersCurrentGroup + 1ULL)) / 2ULL);
			imputedRawPtr = &(imputedRaw[0]);
			//column
			for(std::size_t i = 0; i < nMarkersCurrentGroup; i++)
			{
				//row
				for(std::size_t j = 0; j <= i; j++)
				{
					std::size_t row = (std::size_t)markersThisGroup[j], column = (std::size_t)markersThisGroup[i];
					if(row > column) std::swap(row, column);
					imputedRaw[(i * (i + 1ULL))/2ULL + j] = thetaRawData[(column * (column + 1ULL))/2ULL + row];
				}
			}

			std::string error;
			std::function<void(unsigned long,unsigned long)> imputationProgressFunction = [](unsigned long,unsigned long){};
			if(verbose)
			{
				Rcpp::Rcout << "Starting imputation for group " << *currentGroup << std::endl;
				barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
					imputationProgressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
				{
#ifdef CUSTOM_STATIC_RCPP
					setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#else
					setTxtProgressBar(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#endif
				};
			}
			std::vector<unsigned char> thetaRawDataCopied = imputedRaw;
			unsigned char* originalRawPtr = &(thetaRawDataCopied[0]);
			std::vector<std::pair<int, int> > reportedErrors;
			bool imputationResult = impute(originalRawPtr, imputedRawPtr, levels, NULL, NULL, contiguousIndices, imputationProgressFunction, false, reportedErrors);
			if(std::find(imputedRaw.begin(), imputedRaw.end(), 255) != imputedRaw.end())
			{
				throw std::runtime_error("Internal error in order.cpp");
			}
			if(verbose)
			{
				close(barHandle);
			}

			if(!imputationResult)
			{
				if(reportedErrors.size() > 0)
				{
					std::stringstream ss;
					ss << "Error performing imputation for group " << *currentGroup << ": Unable to impute a value for marker " << reportedErrors[0].first << " and marker "<< reportedErrors[1].second;
					throw std::runtime_error(ss.str().c_str());
				}
				else
				{
					throw std::runtime_error("Internal error");
				}
			}
		}
		else
		{
			imputedRawPtr = &(Rcpp::as<Rcpp::RawVector>(Rcpp::as<Rcpp::S4>(imputedTheta(groupCount)).slot("data"))[0]);
		}
		//Unpack the data into a symmetric matrix
		std::vector<Rbyte> distMatrix(nMarkersCurrentGroup*nMarkersCurrentGroup);
		for(std::size_t i = 0; i < nMarkersCurrentGroup; i++)
		{
			for(std::size_t j = 0; j <= i; j++)
			{
				distMatrix[i * nMarkersCurrentGroup + j] = distMatrix[j * nMarkersCurrentGroup + i] = imputedRawPtr[(i*(i+1ULL))/2ULL + j];
			}
		}
	
		
		std::function<bool(unsigned long, unsigned long)> orderingProgressFunction = [](unsigned long,unsigned long){return false;};
		if(verbose)
		{
			//Only output this text if there was an imputation step, or we're ordering multiple groups
			if(!hasImputedTheta || groupsToOrder.size() > 1) Rcpp::Rcout << "Starting to order group " << *currentGroup << std::endl;
			barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
			orderingProgressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
			{
#ifdef CUSTOM_STATIC_RCPP
				setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#else
				setTxtProgressBar(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#endif
				return false;
			};
		}
		arsaRaw::arsaRawArgs args(levels, currentGroupPermutation);
		args.n = nMarkersCurrentGroup;
		args.rawDist = &(distMatrix[0]);
		args.cool = cool;
		args.temperatureMin = temperatureMin;
		args.nReps = nReps;
		args.progressFunction = orderingProgressFunction;
		args.randomStart = randomStart;
		args.maxMove = maxMove;
		args.effortMultiplier = effortMultiplier;
#ifdef USE_OPENMP
		if(omp_get_max_threads() > 1)
		{
			arsaRaw::arsaRawParallel(args);
		}
		else
#endif
		{
			arsaRaw::arsaRaw(args);
		}

		if(verbose)
		{
			close(barHandle);
		}
		for(std::size_t i = 0; i < nMarkersCurrentGroup; i++) permutation.push_back(markersThisGroup[currentGroupPermutation[i]]+1);
	}
	return Rcpp::wrap(permutation);
END_RCPP
}
