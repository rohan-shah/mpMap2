#ifndef GET_ALLOWABLE_MARKER_PATTERNS_HEADER_GUARD
#define GET_ALLOWABLE_MARKER_PATTERNS_HEADER_GUARD
#include <vector>
#include <map>
#include "unitTypes.hpp"
void getAllowableMarkerPatterns(std::vector<bool>& allowableMarkerPatterns, std::vector<bool>& allowableMarkerPatternsIRIP, std::map<markerEncoding, markerPatternID>& markerPatterns, unsigned int nFounders);
#endif