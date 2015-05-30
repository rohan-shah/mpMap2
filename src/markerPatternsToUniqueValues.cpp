#include "markerPatternsToUniqueValues.h"
void markerPatternsToUniqueValues(markerPatternsToUniqueValuesArgs& args)
{
	for(long markerCounter = 0; markerCounter < args.nMarkers; markerCounter++)
	{
		Rcpp::IntegerMatrix currentMarkerHetData = args.recodedHetData(markerCounter);
		//Set up struct containing data about this marker
		markerData currentPattern(currentMarkerHetData.nrow(), currentMarkerHetData.ncol(), args.nFounders);
		for(long founderCounter = 0; founderCounter < args.nFounders; founderCounter++)
		{
			currentPattern.founderAlleles[founderCounter] = args.recodedFounders(founderCounter, markerCounter);
		}
		for(int i = 0; i < currentMarkerHetData.nrow(); i++)
		{
			for(int j = 0; j < currentMarkerHetData.ncol(); j++)
			{
				currentPattern.hetData(i, j) = currentMarkerHetData(i, j);
			}
		}
		currentPattern.nObservedValues = *std::max_element(&(currentPattern.founderAlleles[0]), &(currentPattern.founderAlleles[args.nFounders]));
		currentPattern.nObservedValues = std::max(currentPattern.nObservedValues, *std::max_element(&(currentPattern.hetData(0,0)), &(currentPattern.hetData(0,0))+currentPattern.hetData.getNRows() * currentPattern.hetData.getNColumns()));
		currentPattern.nObservedValues++;

		currentPattern.computeHash();
		//Have we already seen a marker like this?
		std::map<markerData, markerPatternID>::iterator existingPattern = args.markerPatterns.find(currentPattern);
		//If not, add it to the list of all possible marker patterns, and give it a new ID
		if(existingPattern == args.markerPatterns.end())
		{
			args.markerPatternIDs.push_back(args.markerPatterns.size());
			args.markerPatterns.insert(std::make_pair(currentPattern, args.markerPatterns.size()));
			args.allMarkerPatterns.push_back(currentPattern);
		}
		else
		{
			//If we have, mark this marker as having pattern given by the existing ID
			args.markerPatternIDs.push_back(existingPattern->second);
		}
	}
}