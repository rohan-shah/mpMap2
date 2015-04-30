#include "intercrossingGenerations.h"
bool getIntercrossingGenerations(Rcpp::S4 pedigree, Rcpp::IntegerMatrix finals, int nFounders, std::vector<int>& output)
{
	Rcpp::CharacterVector pedigreeLineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));
	int nPedigreeRows = pedigreeLineNames.length();

	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother")), father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[0]);
	int nFinals = finals.nrow();

	output.clear();
	output.resize(nFinals);
	for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
	{
		Rcpp::CharacterVector::iterator findLineNameInPedigree = std::find(pedigreeLineNames.begin(), pedigreeLineNames.end(), finalNames[finalCounter]);
		if(findLineNameInPedigree == pedigreeLineNames.end()) return false;
		int currentPedRow = std::distance(pedigreeLineNames.begin(), findLineNameInPedigree);
		//First we go backwards through the selfing generations
		while(mother(currentPedRow) == father(currentPedRow))
		{
			currentPedRow = mother(currentPedRow)-1;
			if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;
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
		output[finalCounter] = nAIC - (int)((log(nFounders) / log(2)) + 0.5);
	}
	return true;
}
