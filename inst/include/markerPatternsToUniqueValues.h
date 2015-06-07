#ifndef MARKER_PATTERNS_TO_UNIQUE_VALUES
#define MARKER_PATTERNS_TO_UNIQUE_VALUES
#include <map>
#include <vector>
#include <Rcpp.h>
#include "unitTypes.hpp"
#include "matrices.hpp"
#include "crc32.h"
struct markerData
{
public:
	markerData(int nFounders)
		: hetData(nFounders, nFounders, -1), hashed(false), hash(0)
	{}
	rowMajorMatrix<int> hetData;
	bool operator<(const markerData& other) const
	{
		if(!this->hashed) computeHash();
		if(!other.hashed) other.computeHash();
		if(this->hash < other.hash) return true;
		if(this->hash > other.hash) return false;
		for(int i = 0; i < hetData.getNRows(); i++)
		{
			for(int j = 0; j < hetData.getNColumns(); j++)
			{
				if(this->hetData(i, j) < other.hetData(i, j)) return true;
				if(this->hetData(i, j) > other.hetData(i, j)) return false;
			}
		}
		return false;
	}
	void computeHash() const
	{
		if(!hashed)
		{
			hash = crc32(&(hetData(0,0)), sizeof(int)*hetData.getNRows() * hetData.getNColumns(), hash);
			hashed = true;
		}
	}
	int nObservedValues;
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
	void swap(markerPatternsToUniqueValuesArgs& other)
	{
		markerPatterns.swap(other.markerPatterns);
		markerPatternIDs.swap(other.markerPatternIDs);
		allMarkerPatterns.swap(other.allMarkerPatterns);
		std::swap(nFounders, other.nFounders);
		std::swap(nMarkers, other.nMarkers);
		recodedFounders = other.recodedFounders;
		recodedHetData = other.recodedHetData;
	}
};
void markerPatternsToUniqueValues(markerPatternsToUniqueValuesArgs& args);
#endif