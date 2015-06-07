#include "markerPatternsToUniqueValues.h"
void markerPatternsToUniqueValues(markerPatternsToUniqueValuesArgs& args)
{
	for(long markerCounter = 0; markerCounter < args.nMarkers; markerCounter++)
	{
		Rcpp::IntegerMatrix currentMarkerHetData = args.recodedHetData(markerCounter);
		//Set up struct containing data about this marker
		markerData currentPattern(args.nFounders);
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
					}
				}
			}
		}
		currentPattern.nObservedValues = *std::max_element(&(currentPattern.hetData(0,0)), &(currentPattern.hetData(0,0))+currentPattern.hetData.getNRows() * currentPattern.hetData.getNColumns());
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