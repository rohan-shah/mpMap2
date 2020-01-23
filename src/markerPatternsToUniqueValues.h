#ifndef MARKER_PATTERNS_TO_UNIQUE_VALUES
#define MARKER_PATTERNS_TO_UNIQUE_VALUES
#include <map>
#include <vector>
#include <Rcpp.h>
#include "unitTypes.h"
#include "matrices.h"
#include "crc32.h"
struct markerData
{
public:
	markerData(int nFounders);
	bool operator<(const markerData& other) const;
	void computeHash() const;
	rowMajorMatrix<int> hetData;
	int nObservedValues;
	markerData(const markerData& other)
		:nObservedValues(other.nObservedValues), hashed(other.hashed), hash(other.hash) 
	{
		hetData = other.hetData.copy();
	}
private:
	mutable bool hashed;
	mutable int hash;
};
struct markerPatternsToUniqueValuesArgs
{
public:
	std::map<markerData, markerPatternID> markerPatterns;
	std::vector<markerPatternID> markerPatternIDs;
	std::vector<markerData> allMarkerPatterns;
	int nFounders;
	int nMarkers;
	Rcpp::IntegerMatrix recodedFounders;
	Rcpp::List recodedHetData;
	void swap(markerPatternsToUniqueValuesArgs& other);
};
void markerPatternsToUniqueValues(markerPatternsToUniqueValuesArgs& args);
#endif
