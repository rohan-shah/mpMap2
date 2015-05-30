#include "intercrossingAndSelfingGenerations.h"
bool getIntercrossingAndSelfingGenerations(Rcpp::S4 pedigree, Rcpp::IntegerMatrix finals, int nFounders, std::vector<int>& intercrossing, std::vector<int>& selfing)
{
	Rcpp::CharacterVector pedigreeLineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));
	int nPedigreeRows = pedigreeLineNames.length();

	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother")), father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[0]);
	int nFinals = finals.nrow();

	intercrossing.clear();
	intercrossing.resize(nFinals);

	selfing.clear();
	selfing.resize(nFinals);
	int log2Founders = (int)((log(nFounders) / log(2)) + 0.5);
	for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
	{
		Rcpp::CharacterVector::iterator findLineNameInPedigree = std::find(pedigreeLineNames.begin(), pedigreeLineNames.end(), finalNames[finalCounter]);
		if(findLineNameInPedigree == pedigreeLineNames.end()) return false;
		int currentPedRow = std::distance(pedigreeLineNames.begin(), findLineNameInPedigree);
		int currentSelfing = 0;
		//First we go backwards through the selfing generations
		while(mother(currentPedRow) == father(currentPedRow))
		{
			currentPedRow = mother(currentPedRow)-1;
			if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;
			currentSelfing++;
		}
		//Now the AIC generations
		int nAIC = 0;
		while(mother(currentPedRow) > 0)
		{
			currentPedRow = mother(currentPedRow)-1;
			if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;
			nAIC++;
		}
		//But we've counted the initial log2(nFounders) generations, so subtract that
		intercrossing[finalCounter] = nAIC - log2Founders;
		selfing[finalCounter] = currentSelfing;
	}
	return true;
}
