#include "estimateRfCheckFunnels.h"
#include "estimateRf.h"
#include "getFunnel.h"
#include "orderFunnel.h"
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
void estimateRfCheckFunnels(Rcpp::IntegerMatrix finals, Rcpp::IntegerMatrix founders, Rcpp::S4 pedigree, std::vector<int>& intercrossingGenerations, std::vector<std::string>& warnings, std::vector<std::string>& errors, std::vector<funnelType>& allFunnels)
{
	Rcpp::CharacterVector pedigreeLineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));
	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother"));
	Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[0]);
	Rcpp::CharacterVector markerNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[1]);
	int nFinals = finals.nrow(), nFounders = founders.nrow(), nMarkers = finals.ncol();

	std::vector<long> individualsToCheckFunnels;
	for(long finalCounter = 0; finalCounter < nFinals; finalCounter++)
	{
		individualsToCheckFunnels.clear();

		Rcpp::CharacterVector::iterator findLineName = std::find(pedigreeLineNames.begin(), pedigreeLineNames.end(), finalNames[finalCounter]);
		if(findLineName == pedigreeLineNames.end())
		{
			std::stringstream ss;
			ss << "Unable to find line number " << finalCounter << " named " << finalNames[finalCounter] << " in pedigree";
			throw std::runtime_error(ss.str().c_str());
		}
		int pedigreeRow = (int)std::distance(pedigreeLineNames.begin(), findLineName);
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
				ss << "Attempting to trace pedigree for line " << finalNames[finalCounter] << ": Unable to get funnel for line " << pedigreeLineNames[*i];
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
				ss << "Funnel for line " << pedigreeLineNames[*i] << " contained founders {" << funnel.val[0];
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

			for(std::vector<int>::iterator founderIterator = representedFounders.begin(); founderIterator != representedFounders.end(); founderIterator++)
			{
				//Note that founderIterator comes from representedFounders, which comes from getFunnel - Which returns values starting at 1, not 0. So we have to subtract one. 
				if(finals(finalCounter, markerCounter) == founders((*founderIterator)-1, markerCounter))
				{
					okMarker = true;
					break;
				}
			}
			if(!okMarker)
			{
				std::stringstream ss;
				ss << "Error: Data for marker " << markerNames[markerCounter] << " is impossible for individual " << finalNames[finalCounter] << " with given pedigree\n";
				errors.push_back(ss.str());
			}
		}
		//In this case individualsToCheckFunnels contains one element => getFunnel was only called once => we can reuse the funnel variable
		if(intercrossingGenerations[finalCounter] == 0)
		{
			orderFunnel(&(funnel.val[0]), nFounders);
			allFunnels.push_back(funnel);
		}
	}
}