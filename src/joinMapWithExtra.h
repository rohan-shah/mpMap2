#ifndef JOIN_MAP_WITH_EXTRA_HEADER_GUARD
#define JOIN_MAP_WITH_EXTRA_HEADER_GUARD
#include <Rcpp.h>
struct positionData
{
	Rcpp::List makeUnifiedMap();
	struct chromosomeDescriptor
	{
		int start, end;
		std::string name;
	};
	std::vector<std::string> names;
	std::vector<double> positions;
	std::vector<int> markerIndices;
	std::vector<chromosomeDescriptor> chromosomes;
};
void joinMapWithExtra(Rcpp::List map, Rcpp::List extraPositions, positionData& allPositions);
#endif
