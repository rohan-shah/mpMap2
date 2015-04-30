#include "markerPatternsToUniqueValues.h"
void markerPatternsToUniqueValues(std::map<markerEncoding, markerPatternID>& markerPatterns, std::vector<markerPatternID>& markerPatternIDs, std::vector<markerEncoding>&markerEncodings, int nFounders, int nMarkers, Rcpp::IntegerMatrix& recodedFounders)
{
	for(long markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		int encodedMarkerPattern = 0;
		for(long founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encodedMarkerPattern += (recodedFounders(founderCounter, markerCounter) << 3*founderCounter);
		}
		if(markerPatterns.find(encodedMarkerPattern) == markerPatterns.end())
		{
			markerPatternIDs.push_back(markerPatterns.size());
			markerPatterns.insert(std::make_pair(encodedMarkerPattern, markerPatterns.size()));
			markerEncodings.push_back(encodedMarkerPattern);
		}
		else
		{
			markerPatternIDs.push_back(markerPatterns.find(encodedMarkerPattern)->second);
		}
	}
}