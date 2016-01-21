#include "estimateRFCheckFunnels.h"
#include "estimateRF.h"
#include "getFunnel.h"
#include "orderFunnel.h"
#include "matrices.hpp"
#include "sortPedigreeLineNames.h"
void getAICParentLines(Rcpp::IntegerVector& mother, Rcpp::IntegerVector& father, long pedigreeRow, int intercrossingGenerations, std::vector<long>& individualsToCheckFunnels)
{
	//The lines that we currently need to check goes in individualsToCheckFunnels
	individualsToCheckFunnels.clear();
	//convert initial ID into a row in the pedigree
	int currentPedigreeRow = pedigreeRow;
	//step back through all the selfing generations, if they're explicitly listed
	while(mother(currentPedigreeRow) == father(currentPedigreeRow))
	{
		currentPedigreeRow = mother(currentPedigreeRow)-1;
	}
	individualsToCheckFunnels.push_back(currentPedigreeRow);

	//step back through the intercrossing generations
	std::vector<long> nextGenerationToCheck;
	for(;intercrossingGenerations > 0; intercrossingGenerations--)
	{
		nextGenerationToCheck.clear();
		//For everything that currently needs to be checked, go up one generation. 
		for(std::vector<long>::iterator i = individualsToCheckFunnels.begin(); i != individualsToCheckFunnels.end(); i++)
		{
			nextGenerationToCheck.push_back(mother(*i)-1);
			nextGenerationToCheck.push_back(father(*i)-1);
		}
		//swap vectors
		individualsToCheckFunnels.swap(nextGenerationToCheck);
	}
}
/* This function specifically checks whether the observed data is consistent with the *pedigree*. It assumes that every observed value in the finals is already valid - That is, every observed value contained in the finals is also listed as a possibility in the hetData object
*/
void estimateRFCheckFunnels(Rcpp::IntegerMatrix finals, Rcpp::IntegerMatrix founders, Rcpp::List hetData, Rcpp::S4 pedigree, std::vector<int>& intercrossingGenerations, std::vector<std::string>& warnings, std::vector<std::string>& errors, std::vector<funnelType>& allFunnels, std::vector<funnelType>& lineFunnels)
{
	Rcpp::CharacterVector pedigreeLineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));

	//We make a copy of the pedigree line names and sort it (otherwise the std::find relating to pedigreeLineNames is prohibitive)
	std::vector<pedigreeLineStruct> sortedLineNames;
	sortPedigreeLineNames(pedigreeLineNames, sortedLineNames);

	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother"));
	Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[0]);
	Rcpp::CharacterVector markerNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[1]);
	int nFinals = finals.nrow(), nFounders = founders.nrow(), nMarkers = finals.ncol();

	xMajorMatrix<int> foundersToMarkerAlleles(nFounders, nFounders, nMarkers, -1);
	for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		Rcpp::IntegerMatrix currentMarkerHetData = hetData(markerCounter);
		for(int founderCounter1 = 0; founderCounter1 < nFounders; founderCounter1++)
		{
			for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
			{
				int markerAllele1 = founders(founderCounter1, markerCounter);
				int markerAllele2 = founders(founderCounter2, markerCounter);
				for(int hetDataRowCounter = 0; hetDataRowCounter < currentMarkerHetData.nrow(); hetDataRowCounter++)
				{
					if(markerAllele1 == currentMarkerHetData(hetDataRowCounter, 0) && markerAllele2 == currentMarkerHetData(hetDataRowCounter, 1))
					{
						foundersToMarkerAlleles(founderCounter1, founderCounter2, markerCounter) = currentMarkerHetData(hetDataRowCounter, 2);
					}
				}
			}
		}
	}
	std::vector<long> individualsToCheckFunnels;
	for(long finalCounter = 0; finalCounter < nFinals; finalCounter++)
	{
		individualsToCheckFunnels.clear();
		std::string currentLineName = Rcpp::as<std::string>(finalNames(finalCounter));

		std::vector<pedigreeLineStruct>::iterator findLineName = std::lower_bound(sortedLineNames.begin(), sortedLineNames.end(), pedigreeLineStruct(currentLineName, -1));
		if(findLineName == sortedLineNames.end() || findLineName->lineName != currentLineName)
		{
			std::stringstream ss;
			ss << "Unable to find line number " << finalCounter << " named " << finalNames(finalCounter) << " in pedigree";
			throw std::runtime_error(ss.str().c_str());
		}
		int pedigreeRow = findLineName->index;
		if(intercrossingGenerations[finalCounter] == 0)
		{
			individualsToCheckFunnels.push_back(pedigreeRow);
		}
		else
		{
			getAICParentLines(mother, father, pedigreeRow, intercrossingGenerations[finalCounter], individualsToCheckFunnels);
		}
		//Now we know the lines for which we need to check the funnels from the pedigree (note: We don't necessarily have genotype data for all of these, it's purely a pedigree check)
		//This vector lists all the founders that are ancestors of the current line. This may comprise any number - E.g. if we have an AIC line descended from funnels 1,2,1,2 and 2,3,2,3 then this vector is going it contain 1,2,3
		std::vector<int> representedFounders;
		//Fixed length arrays to store funnels. If we have less than 16 founders then part of this is garbage and we don't use that bit....
		funnelType funnel, copiedFunnel;
		for(std::vector<long>::iterator i = individualsToCheckFunnels.begin(); i != individualsToCheckFunnels.end(); i++)
		{
			try
			{
				getFunnel(*i, mother, father, &(funnel.val[0]), nFounders);
			}
			catch(...)
			{
				std::stringstream ss;
				ss << "Attempting to trace pedigree for line " << finalNames(finalCounter) << ": Unable to get funnel for line " << pedigreeLineNames(*i);
				errors.push_back(ss.str());
				continue;
			}
			//insert these founders into the vector containing all the represented founders
			representedFounders.insert(representedFounders.end(), &(funnel.val[0]), &(funnel.val[0]) + nFounders);
			//Copy the funnel 
			memcpy(&copiedFunnel, &funnel, sizeof(funnelType));
			std::sort(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders);
			if(std::unique(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders) != &(copiedFunnel.val[0]) + nFounders)
			{
				std::stringstream ss;
				ss << "Funnel for line " << pedigreeLineNames(*i) << " contained founders {" << funnel.val[0];
				warnings.push_back(ss.str());
				if(nFounders == 2)
				{
					ss << ", " << funnel.val[1] << "}";
				}
				else if(nFounders == 4)
				{
					ss << ", " << funnel.val[1] << ", " << funnel.val[2] << ", " << funnel.val[3] << "}";
				}
				else
				{
					ss << ", " << funnel.val[1] << ", " << funnel.val[2] << ", " << funnel.val[3] << ", " << funnel.val[4] << ", " << funnel.val[5] << ", " << funnel.val[6] << ", " << funnel.val[7]<< "}";
				}
				ss << ". Did you intend to use all " << nFounders << " founders?";
			}
			else
			{
				allFunnels.push_back(funnel);
			}
		}
		//remove duplicates in representedFounders
		std::sort(representedFounders.begin(), representedFounders.end());
		representedFounders.erase(std::unique(representedFounders.begin(), representedFounders.end()), representedFounders.end());
		//Not having all the founders in the input funnels is more serious if it causes the observed marker data to be impossible. So check for this.
		for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
		{
			bool okMarker = false;
			//If observed value is an NA then than's ok, continue
			int value = finals(finalCounter, markerCounter);
			if(value == NA_INTEGER) continue;

			for(std::vector<int>::iterator founderIterator1 = representedFounders.begin(); founderIterator1 != representedFounders.end(); founderIterator1++)
			{
				for(std::vector<int>::iterator founderIterator2 = representedFounders.begin(); founderIterator2 != representedFounders.end(); founderIterator2++)
				{
					//Note that founderIterator comes from representedFounders, which comes from getFunnel - Which returns values starting at 1, not 0. So we have to subtract one. 
					if(finals(finalCounter, markerCounter) == foundersToMarkerAlleles((*founderIterator1)-1, (*founderIterator2)-1, markerCounter))
					{
						okMarker = true;
						break;
					}
				}
			}
			if(!okMarker)
			{
				std::stringstream ss;
				ss << "Error: Data for marker " << markerNames(markerCounter) << " is impossible for individual " << finalNames(finalCounter) << " with given pedigree\n";
				errors.push_back(ss.str());
				if(errors.size() > 1000) return;
			}
		}
		//In this case individualsToCheckFunnels contains one element => getFunnel was only called once => we can reuse the funnel variable
		if(intercrossingGenerations[finalCounter] == 0)
		{
			orderFunnel(&(funnel.val[0]), nFounders);
			lineFunnels.push_back(funnel);
		}
		else
		{
			//Add a dummy value in lineFunnel
			for(int i = 0; i < 8; i++) funnel.val[i] = 0;
			lineFunnels.push_back(funnel);
		}
	}
}
