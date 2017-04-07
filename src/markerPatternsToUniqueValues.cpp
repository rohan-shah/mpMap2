#include "markerPatternsToUniqueValues.h"
markerData::markerData(int nFounders)
	: hetData(nFounders, nFounders, NA_INTEGER), hashed(false), hash(0)
{}

bool markerData::operator<(const markerData& other) const
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
void markerData::computeHash() const
{
	if(!hashed)
	{
		hash = crc32(&(hetData(0,0)), sizeof(int)*hetData.getNRows() * hetData.getNColumns(), hash);
		hashed = true;
	}
}
void markerPatternsToUniqueValues(markerPatternsToUniqueValuesArgs& args)
{
	for(long markerCounter = 0; markerCounter < args.nMarkers; markerCounter++)
	{
		Rcpp::IntegerMatrix currentMarkerHetData = args.recodedHetData(markerCounter);
		//Set up struct containing data about this marker
		markerData currentPattern(args.nFounders);
		//Variable to check that at least one hetData entry is set in the loop. This could fail in the case that the hetData struct has zero rows, which is allowed in the case of a marker with NA founders. 
		bool hasHetdataEntry = false;
		for(int founderCounter1 = 0; founderCounter1 < args.nFounders; founderCounter1++)
		{
			for(int founderCounter2 = 0; founderCounter2 < args.nFounders; founderCounter2++)
			{
				int markerAllele1 = args.recodedFounders(founderCounter1, markerCounter);
				int markerAllele2 = args.recodedFounders(founderCounter2, markerCounter);
				for(int hetDataRowCounter = 0; hetDataRowCounter < currentMarkerHetData.nrow(); hetDataRowCounter++)
				{
					if(currentMarkerHetData(hetDataRowCounter, 0) == markerAllele1 && currentMarkerHetData(hetDataRowCounter, 1) == markerAllele2)
					{
						currentPattern.hetData(founderCounter1, founderCounter2) = currentMarkerHetData(hetDataRowCounter, 2);
						hasHetdataEntry = true;
					}
				}
			}
		}
		if(!hasHetdataEntry)
		{
			currentPattern.nObservedValues = 0;
		}
		else
		{
			currentPattern.nObservedValues = *std::max_element(&(currentPattern.hetData(0,0)), &(currentPattern.hetData(0,0))+currentPattern.hetData.getNRows() * currentPattern.hetData.getNColumns());
			currentPattern.nObservedValues++;
		}

		currentPattern.computeHash();
		//Have we already seen a marker like this?
		std::map<markerData, markerPatternID>::iterator existingPattern = args.markerPatterns.find(currentPattern);
		//If not, add it to the list of all possible marker patterns, and give it a new ID
		if(existingPattern == args.markerPatterns.end())
		{
			args.markerPatternIDs.push_back((int)args.markerPatterns.size());
			args.markerPatterns.insert(std::make_pair(currentPattern, args.markerPatterns.size()));
			args.allMarkerPatterns.emplace_back(std::move(currentPattern));
		}
		else
		{
			//If we have, mark this marker as having pattern given by the existing ID
			args.markerPatternIDs.push_back(existingPattern->second);
		}
	}
}
void markerPatternsToUniqueValuesArgs::swap(markerPatternsToUniqueValuesArgs& other)
{
	markerPatterns.swap(other.markerPatterns);
	markerPatternIDs.swap(other.markerPatternIDs);
	allMarkerPatterns.swap(other.allMarkerPatterns);
	std::swap(nFounders, other.nFounders);
	std::swap(nMarkers, other.nMarkers);
	recodedFounders = other.recodedFounders;
	recodedHetData = other.recodedHetData;
}

