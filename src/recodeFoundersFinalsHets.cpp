#include "recodeFoundersFinalsHets.h"
/*
	Re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles at that particular marker. The maximum number of alleles across all the markers is recorded and output. 
*/
void recodeFoundersFinalsHets(recodeDataStruct& inputs)
{
	inputs.maxAlleles = 0;
	long nMarkers = inputs.finals.ncol();
	long nFounders = inputs.founders.nrow(), nFinals = inputs.finals.nrow();
	
	std::vector<int> finalAllowedValues;
	std::map<int, int> founderTranslations, finalTranslations;
	for(long markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		founderTranslations.clear();
		finalAllowedValues.clear();	
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
				inputs.recodedFounders(founderCounter, markerCounter) = founderTranslations.size();
				founderTranslations.insert(std::make_pair(oldValue, founderTranslations.size()));
			}
			else
			{
				inputs.recodedFounders(founderCounter, markerCounter) = lookup->second;
			}
		}
		inputs.maxAlleles = std::max(inputs.maxAlleles, (unsigned int)founderTranslations.size());
		//now translate hetData to sequential values
		//Start off by translating the first two columns. Then gather up all the values in the third column that don't correspond to a homozygote. 
		//But if they DO correspond to a homozygote, translate immediately
		for(int i = 0; i < currentMarkerHetData.nrow(); i++)
		{
			for(int j = 0; j < 2; j++)
			{
				recodedCurrentMarkerHetData(i, j) = founderTranslations.find(currentMarkerHetData(i, j))->second;
			}
			if(currentMarkerHetData(i, 0) != currentMarkerHetData(i, 1)) finalAllowedValues.push_back(currentMarkerHetData(i, 2));
			else
			{
				recodedCurrentMarkerHetData(i, 2) = recodedCurrentMarkerHetData(i, 0);
			}
		}
		//The homozygotes get translated as-is, the hetrozygotes are given sequential labels afterwards. 
		finalAllowedValues.erase(std::unique(finalAllowedValues.begin(), finalAllowedValues.end()), finalAllowedValues.end());
		for(std::size_t allowedValueCounter = 0; allowedValueCounter < finalAllowedValues.size(); allowedValueCounter++)
		{
			finalTranslations.insert(std::make_pair(finalAllowedValues[allowedValueCounter], (int)allowedValueCounter+nFounders));
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
				inputs.recodedFinals(finalCounter, markerCounter) = founderTranslations.find(inputs.finals(finalCounter, markerCounter))->second;
			}
			else inputs.recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
		}
		inputs.recodedHetData(markerCounter) = recodedCurrentMarkerHetData;
	}
}
