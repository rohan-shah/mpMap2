#include "impute.h"
#include <vector>
#include <set>
#include <math.h>
#include <limits>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
struct rowColumnDifference
{
public:
	rowColumnDifference(bool isRow, int index, double difference)
		: isRow(isRow), index(index), difference(difference)
	{}
	bool isRow;
	int index;
	double difference;
};
struct differenceSetComp
{
	bool operator()(const rowColumnDifference& a, const rowColumnDifference& b) const
	{
		return a.difference < b.difference;
	}
};
template<bool hasLOD, bool hasLKHD> bool imputeInternal(const unsigned char* originalTheta, unsigned char* imputedTheta, std::vector<double>& levels, double* lod, double* lkhd, std::vector<int>& markersThisGroup, std::function<void(unsigned long, unsigned long)> statusFunction, bool allErrors, std::vector<std::pair<int, int> >& reportedErrors)
{
	unsigned long done = 0;
	unsigned long total = (unsigned long)markersThisGroup.size();
	bool hasError = false;
	int nLevels = (int)levels.size(), nMarkersThisGroup = (int)markersThisGroup.size();
	std::vector<double> absoluteDifferences(0x100*0x100, 0);
	for(int i = 0; i < nLevels; i++)
	{
		for(int j = 0; j < nLevels; j++)
		{
			absoluteDifferences[i * 0x100 + j] = fabs(levels[i] - levels[j]);
		}
	}
	std::vector<int> copiedMarkersThisGroup = markersThisGroup;
	std::sort(copiedMarkersThisGroup.begin(), copiedMarkersThisGroup.end());
	std::vector<double> similarityMatrix(nMarkersThisGroup * nMarkersThisGroup, 0);
	std::vector<int> usablePoints(nMarkersThisGroup * nMarkersThisGroup, 0);
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		std::vector<int> table(0x100 * 0x100, 0);
#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
#endif
		//First marker
		for(unsigned long long i = 0; i < (unsigned long long)nMarkersThisGroup; i++)
		{
			//Second marker
			for(unsigned long long j = 0; j < i; j++)
			{
				std::fill(table.begin(), table.end(), 0);
				unsigned long long usableSum = 0;
				for(unsigned long long k = 0; k <= j; k++)
				{
					table[originalTheta[i * (i + 1ULL) / 2ULL + k] * 0x100 + originalTheta[j * (j + 1ULL) / 2ULL + k]]++;
					usableSum += (originalTheta[i * (i + 1ULL) / 2ULL + k] != 0xff && originalTheta[j * (j + 1ULL) / 2ULL + k] != 0xff);
				}
				for(unsigned long long k = j+i; k <= i; k++)
				{
					table[originalTheta[i * (i + 1ULL) / 2ULL + k] * 0x100 + originalTheta[k * (k + 1ULL) / 2ULL + j]]++;
					usableSum += (originalTheta[i * (i + 1ULL) / 2ULL + k] != 0xff && originalTheta[k * (k + 1ULL) / 2ULL + j] != 0xff);
				}
				for(unsigned long long k = i+1; k < (unsigned long long)nMarkersThisGroup; k++)
				{
					table[originalTheta[k * (k + 1ULL) / 2ULL + i] * 0x100 + originalTheta[k * (k + 1ULL) / 2ULL + j]]++;
					usableSum += (originalTheta[k * (k + 1ULL) / 2ULL + i] != 0xff && originalTheta[k * (k + 1ULL) / 2ULL + j] != 0xff);
				}
				double sum = 0;
				for(int k1 = 0; k1 < nLevels; k1++)
				{
					for(int k2 = 0; k2 < nLevels; k2++)
					{
						sum += table[k1 * 0x100 + k2] * absoluteDifferences[k1 * 0x100 + k2];
					}
				}
				similarityMatrix[j * nMarkersThisGroup + i] = similarityMatrix[i * nMarkersThisGroup + j] = sum;
				usablePoints[i*nMarkersThisGroup + j] = usablePoints[j * nMarkersThisGroup + i] = usableSum;
			}
		}
	}
	typedef std::set<rowColumnDifference, differenceSetComp> differenceSetType;
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		differenceSetType differenceSet;
#ifdef _OPENMP
		#pragma omp for schedule(dynamic)
#endif
		//row
		for(unsigned long long i = 0; i < (unsigned long long)nMarkersThisGroup; i++)
		{
			int marker1 = copiedMarkersThisGroup[i];
			//column
			for(unsigned long long j = i; j < nMarkersThisGroup; j++)
			{
				int marker2 = copiedMarkersThisGroup[j];
				unsigned long long copiedMarker1 = marker1;
				unsigned long long copiedMarker2 = marker2;
				if(copiedMarker1 > copiedMarker2) std::swap(copiedMarker1, copiedMarker2);
				//record whether there's a missing value for marker1 anywhere.
				float val = originalTheta[(copiedMarker2 * (copiedMarker2 + 1ULL))/2ULL + copiedMarker1];
				unsigned char& toReplace = imputedTheta[(copiedMarker2 * (copiedMarker2 + 1ULL))/2ULL + copiedMarker1];
				if(val == 0xff && toReplace == 0xff)
				{
					differenceSet.clear();
					//There's a missing value for (i, j). So find the closest matching column or column, and overwrite with the value from that column.
					//i is the row, j is the column, k is the candidate other marker to use for the imputation
					for(unsigned long long k = 0; k < nMarkersThisGroup; k++)
					{
						if(k == i || k == j) continue;
						int marker3 = copiedMarkersThisGroup[k];
						if(usablePoints[i * nMarkersThisGroup + k]) differenceSet.insert(rowColumnDifference(true, k, similarityMatrix[i * nMarkersThisGroup + k]));
						if(usablePoints[j * nMarkersThisGroup + k]) differenceSet.insert(rowColumnDifference(false, k, similarityMatrix[j * nMarkersThisGroup + k]));
					}
					//go through the other markers from most similar to least similar, looking for something which has a value here. So marker3 is the marker which is similar to marker1
					bool replacementFound = false;
					for(differenceSetType::iterator marker3 = differenceSet.begin(); marker3 != differenceSet.end(); marker3++)
					{
						unsigned long long pair2Row, pair2Column;
						if(marker3->isRow)
						{
							pair2Row = marker3->index;
							pair2Column = marker2;
						}
						else
						{
							pair2Row = marker1;
							pair2Column = marker3->index;
						}
						if(pair2Row > pair2Column) std::swap(pair2Row, pair2Column);
						const unsigned char& newValue = originalTheta[(pair2Column * (pair2Column + 1ULL))/2ULL + pair2Row];
						if(newValue != 0xff)
						{
							toReplace = newValue;
							if(hasLOD) lod[(copiedMarker2 * (copiedMarker2 + 1ULL))/2ULL + copiedMarker1] = lod[(pair2Column * (pair2Column + 1ULL))/2ULL + pair2Row];
							if(hasLKHD) lkhd[(copiedMarker2 * (copiedMarker2 + 1ULL))/2ULL + copiedMarker1] = lkhd[(pair2Column * (pair2Column + 1ULL))/2ULL + pair2Row];
							//We only need to copy the value from the best other similar marker, so we can break here
							replacementFound = true;
							break;
						}
					}
					if(!replacementFound)
					{
						//critical section, due to assignment to the shared error variable. 
#ifdef _OPENMP
						#pragma omp critical
#endif
						{
							reportedErrors.push_back(std::make_pair(marker1+1, marker2+1));
							hasError = true;
#ifndef _OPENMP
							//We can't return early if we're using openmp
							if(!allErrors) return false;
#endif
						}
					}
				}
			}
#ifdef _OPENMP
			#pragma omp critical
#endif
			{
				done++;
			}
#ifdef _OPENMP
			if(omp_get_thread_num() == 0)
#endif
			{
				statusFunction(done, total);
			}
		}
	}
	return !hasError;
}
bool impute(const unsigned char* originalTheta, unsigned char* imputedTheta, std::vector<double>& thetaLevels, double* lod, double* lkhd, std::vector<int>& markers, std::function<void(unsigned long, unsigned long)> statusFunction, bool allErrors, std::vector<std::pair<int, int> >& reportedErrors)
{
	if(lod != NULL && lkhd != NULL)
	{
		return imputeInternal<true, true>(originalTheta, imputedTheta, thetaLevels, lod, lkhd, markers, statusFunction, allErrors, reportedErrors);
	}
	else if(lod != NULL && lkhd == NULL)
	{
		return imputeInternal<true, false>(originalTheta, imputedTheta, thetaLevels, lod, lkhd, markers, statusFunction, allErrors, reportedErrors);
	}
	else if(lod == NULL && lkhd != NULL)
	{
		return imputeInternal<false, true>(originalTheta, imputedTheta, thetaLevels, lod, lkhd, markers, statusFunction, allErrors, reportedErrors);
	}
	else
	{
		return imputeInternal<false, false>(originalTheta, imputedTheta, thetaLevels, lod, lkhd, markers, statusFunction, allErrors, reportedErrors);
	}
}
SEXP imputeGroup(SEXP mpcrossLG_sexp, SEXP verbose_sexp, SEXP group_sexp, SEXP allErrors_sexp)
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

	bool allErrors;
	try
	{
		allErrors = Rcpp::as<bool>(allErrors_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input allErrors must be a boolean");
	}

	Rcpp::List verboseList;
	bool verbose;
	int progressStyle;
	try
	{
		verboseList = Rcpp::as<Rcpp::List>(verbose_sexp);
		verbose = Rcpp::as<bool>(verboseList("verbose"));
		progressStyle = Rcpp::as<int>(verboseList("progressStyle"));
	}
	catch(...)
	{
		throw std::runtime_error("Input verbose must be a boolean or a list with entries verbose and progressStyle");
	}

	Rcpp::RawVector copiedTheta(((unsigned long long)markersCurrentGroup.size()*((unsigned long long)markersCurrentGroup.size() + 1ULL))/2ULL);
	//column
	for(unsigned long long marker1Counter = 0; marker1Counter < markersCurrentGroup.size(); marker1Counter++)
	{
		//row
		for(unsigned long long marker2Counter = 0; marker2Counter <= marker1Counter; marker2Counter++)
		{
			copiedTheta[(marker1Counter*(marker1Counter+1ULL))/2ULL + marker2Counter] = thetaData[((unsigned long long)markersCurrentGroup[marker1Counter] *((unsigned long long)markersCurrentGroup[marker1Counter] + 1ULL))/2ULL + (unsigned long long)markersCurrentGroup[marker2Counter]];
			if(copiedLodPtr) copiedLodPtr[(marker1Counter*(marker1Counter+1))/2 + marker2Counter] = lodS4Data[((unsigned long long)markersCurrentGroup[marker1Counter] *(markersCurrentGroup[marker1Counter] + 1))/2 + markersCurrentGroup[marker2Counter]];
			if(copiedLkhdPtr) copiedLkhdPtr[(marker1Counter*(marker1Counter+1))/2 + marker2Counter] = lkhdS4Data[((unsigned long long)markersCurrentGroup[marker1Counter] *((unsigned long long)markersCurrentGroup[marker1Counter] + 1ULL))/2ULL + (unsigned long long)markersCurrentGroup[marker2Counter]];
		}
	}
	Rcpp::RawVector originalTheta(copiedTheta.size());
	memcpy(&(originalTheta[0]), &(copiedTheta[0]), sizeof(Rbyte) * copiedTheta.size());

	std::function<void(unsigned long, unsigned long)> progressFunction = [](unsigned long, unsigned long){};
	Rcpp::Function txtProgressBar("txtProgressBar"), setTxtProgressBar("setTxtProgressBar"), close("close"), condition = Rcpp::Environment::namespace_env("mpMap2").get("condition"), signalCondition("signalCondition");
	Rcpp::RObject barHandle;
	if(verbose)
	{
		Rcpp::Rcout << "Starting imputation for group " << group << std::endl;
		barHandle = txtProgressBar(Rcpp::Named("style") = progressStyle, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
		progressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
		{
#ifdef CUSTOM_STATIC_RCPP
			setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#else
			setTxtProgressBar(barHandle, (int)((double)(1000*done) / (double)totalSteps));
#endif
		};
	}
	//Now overwrite the markersCurrentGroup vector with consecutive numbers. Because we've extracted a subset of the matrix into its own memory. 
	for(int i = 0; i < (int)markersCurrentGroup.size(); i++)
	{
		markersCurrentGroup[i] = i;
	}
	std::string error;
	unsigned char* originalThetaPtr = &(originalTheta[0]);
	unsigned char* imputedThetaPtr = &(copiedTheta[0]);
	std::vector<std::pair<int, int> > reportedErrors;
	bool ok = impute(originalThetaPtr, imputedThetaPtr, levels, copiedLodPtr, copiedLkhdPtr, markersCurrentGroup, progressFunction, allErrors, reportedErrors);
	if(!ok)
	{
		if(reportedErrors.size() > 0)
		{
			if(allErrors)
			{
				Rcpp::IntegerMatrix convertedImputationErrors(reportedErrors.size(), 2);
				for(std::size_t i = 0; i < (std::size_t)reportedErrors.size(); i++)
				{
					convertedImputationErrors(i, 0) = reportedErrors[i].first;
					convertedImputationErrors(i, 1) = reportedErrors[i].second;
				}
				signalCondition(condition("imputationErrors", convertedImputationErrors));
			}
			std::stringstream ss;
			ss << "Error performing imputation for group " << group << ": Unable to impute a value for marker " << reportedErrors[0].first << " and marker "<< reportedErrors[1].second;
			throw std::runtime_error(ss.str().c_str());
		}
		else
		{
			throw std::runtime_error("Internal error");
		}
	}
	if(verbose)
	{
		close(barHandle);
	}
	return Rcpp::List::create(Rcpp::Named("theta") = copiedTheta, Rcpp::Named("lod") = copiedLod, Rcpp::Named("lkhd") = copiedLkhd);
END_RCPP
}
