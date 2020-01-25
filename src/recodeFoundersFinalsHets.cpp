#include "recodeFoundersFinalsHets.h"
/*
	Re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles at that particular marker. The maximum number of alleles across all the markers is recorded and output. 
*/
void recodeFoundersFinalsHets(recodeDataStruct& inputs)
{
	inputs.maxAlleles = 0;
	long nMarkers = inputs.finals.ncol();
	long nFounders = inputs.founders.nrow(), nFinals = inputs.finals.nrow();
	
	std::vector<int> hetValues;
	std::map<int, int> founderTranslations, finalTranslations;
	if(Rcpp::as<Rcpp::List>(inputs.hetData).size() != nMarkers)
	{
		throw std::runtime_error("Internal error");
	}
	for(long markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		founderTranslations.clear();
		hetValues.clear();	
		finalTranslations.clear();
		
		Rcpp::IntegerMatrix currentMarkerHetData = Rcpp::as<Rcpp::IntegerMatrix>(Rcpp::as<Rcpp::List>(inputs.hetData)(markerCounter));
		Rcpp::IntegerMatrix recodedCurrentMarkerHetData(currentMarkerHetData.nrow(), currentMarkerHetData.ncol());
		//map to hold the translation from old values to recoded values (0 - (n-1))
		std::map<int, int>::iterator lookup;
		for(long founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			int oldValue = inputs.founders(founderCounter, markerCounter);
			lookup = founderTranslations.find(oldValue);
			if(lookup == founderTranslations.end()) 
			{
				int newValue = (int)founderTranslations.size();
				inputs.recodedFounders(founderCounter, markerCounter) = newValue;
				founderTranslations.insert(std::make_pair(oldValue, newValue));
				finalTranslations.insert(std::make_pair(oldValue, newValue));
			}
			else
			{
				inputs.recodedFounders(founderCounter, markerCounter) = lookup->second;
			}
		}
		unsigned int nFounderAlleles = (unsigned int)founderTranslations.size();
		//now translate hetData to sequential values
		//Start off by translating the first two columns. Then gather up all the values in the third column that don't correspond to a homozygote. 
		//But if the first two columns indicate that we're dealing with a het, translate the third column immediately
		for(int i = 0; i < currentMarkerHetData.nrow(); i++)
		{
			for(int j = 0; j < 2; j++)
			{
				std::map<int, int>::iterator searchResult = founderTranslations.find(currentMarkerHetData(i, j));
				if(searchResult != founderTranslations.end()) recodedCurrentMarkerHetData(i, j) = searchResult->second;
				else throw std::invalid_argument("Invalid genotype encoding");
			}
			if(currentMarkerHetData(i, 0) != currentMarkerHetData(i, 1)) hetValues.push_back(currentMarkerHetData(i, 2));
			else
			{
				recodedCurrentMarkerHetData(i, 2) = recodedCurrentMarkerHetData(i, 0);
			}
		}
		//The homozygotes get translated as-is, the heterozygotes are given sequential labels afterwards. 
		std::sort(hetValues.begin(), hetValues.end());
		hetValues.erase(std::unique(hetValues.begin(), hetValues.end()), hetValues.end());
		inputs.maxAlleles = std::max(inputs.maxAlleles, (unsigned int)(nFounderAlleles + hetValues.size()));
		for(std::size_t allowedValueCounter = 0; allowedValueCounter < hetValues.size(); allowedValueCounter++)
		{
			finalTranslations.insert(std::make_pair(hetValues[allowedValueCounter], (int)allowedValueCounter+nFounderAlleles));
		}
		for(int i = 0; i < currentMarkerHetData.nrow(); i++)
		{
			if(recodedCurrentMarkerHetData(i, 0) != recodedCurrentMarkerHetData(i, 1)) recodedCurrentMarkerHetData(i, 2) = finalTranslations.find(currentMarkerHetData(i, 2))->second;
		}
		
		//now translate finals
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if(inputs.finals(finalCounter, markerCounter) != NA_INTEGER)
			{
				std::map<int, int>::iterator i = finalTranslations.find(inputs.finals(finalCounter, markerCounter));
				if(i != finalTranslations.end()) inputs.recodedFinals(finalCounter, markerCounter) = i->second;
				else throw std::invalid_argument("Invalid genotype encoding");
			}
			else inputs.recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
		}
		inputs.recodedHetData(markerCounter) = recodedCurrentMarkerHetData;
	}
}
