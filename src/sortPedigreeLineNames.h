#ifndef SORT_PEDIGREE_LINE_NAMES_HEADER_GUARD
#define SORT_PEDIGREE_LINE_NAMES_HEADER_GUARD
#include <string>
#include <vector>
#include <Rcpp.h>
struct pedigreeLineStruct
{
public:
	pedigreeLineStruct(std::string lineName, int index)
		: lineName(lineName), index(index)
	{}
	bool operator<(const pedigreeLineStruct& other) const
	{
		return lineName < other.lineName;
	}
	std::string lineName;
	int index;
};
void sortPedigreeLineNames(Rcpp::CharacterVector pedigreeLineNames, std::vector<pedigreeLineStruct>& sortedPedigreeLineNames);
#endif
