#ifndef MARKER_PATTERNS_TO_UNIQUE_VALUES
#define MARKER_PATTERNS_TO_UNIQUE_VALUES
#include <map>
#include <vector>
#include <Rcpp.h>
#include "unitTypes.hpp"
void markerPatternsToUniqueValues(std::map<markerEncoding, markerPatternID>& markerPatterns, std::vector<markerPatternID>& markerPatternIDs, std::vector<markerEncoding>&markerEncodings, int nFounders, int nMarkers, Rcpp::IntegerMatrix& recodedFounders);
#endif