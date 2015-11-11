#include "sortPedigreeLineNames.h"
void sortPedigreeLineNames(Rcpp::CharacterVector pedigreeLineNames, std::vector<pedigreeLineStruct>& sortedPedigreeLineNames)
{
	sortedPedigreeLineNames.clear();
	sortedPedigreeLineNames.reserve(pedigreeLineNames.size());

	for(int i = 0; i < pedigreeLineNames.size(); i++)
	{
		sortedPedigreeLineNames.push_back(pedigreeLineStruct(Rcpp::as<std::string>(pedigreeLineNames[i]), i));
	}
	std::sort(sortedPedigreeLineNames.begin(), sortedPedigreeLineNames.end());
}
