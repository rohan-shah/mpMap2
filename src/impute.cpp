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
