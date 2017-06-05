#include "getAllFunnelsIncAIC.h"
#include "intercrossingAndSelfingGenerations.h"
#include "estimateRF.h"
#include "getFunnel.h"
#include "orderFunnel.h"
#include "estimateRFCheckFunnels.h"
SEXP getAllFunnelsIncAIC(SEXP geneticData_sexp, SEXP standardise_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData;
	try
	{
		geneticData = Rcpp::as<Rcpp::S4>(geneticData_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input must be an S4 object of class geneticData");
	}
	
	Rcpp::IntegerMatrix founders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@founders must be an integer matrix");
	}

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@finals must be an integer matrix");
	}
	Rcpp::CharacterVector finalNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(finals.attr("dimnames"))[0]);

	Rcpp::S4 pedigree;
	try
	{
		pedigree = Rcpp::as<Rcpp::S4>(geneticData.slot("pedigree"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@finals must be an integer matrix");
	}
	Rcpp::CharacterVector pedigreeLineNames = Rcpp::as<Rcpp::CharacterVector>(pedigree.slot("lineNames"));

	bool standardise;
	try
	{
		standardise = Rcpp::as<bool>(standardise_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input standardise must be a boolean");
	}

	int nFounders = founders.nrow();
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(pedigree, finals, nFounders, intercrossingGenerations, selfingGenerations);

	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother"));
	Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	std::vector<funnelType> allFunnels;
	funnelType funnel;

	for(int lineCounter = 0; lineCounter < finalNames.size(); lineCounter++)
	{
		Rcpp::CharacterVector::iterator i = std::find(pedigreeLineNames.begin(), pedigreeLineNames.end(), finalNames[lineCounter]);
		if(i == pedigreeLineNames.end()) throw std::runtime_error("Internal error");
		std::size_t pedigreeRow = std::distance(pedigreeLineNames.begin(), i);
		if(intercrossingGenerations[lineCounter] == 0)
		{
			if(i == pedigreeLineNames.end())
			{
				std::stringstream ss;
				ss << "Unable to find line named " << finalNames[lineCounter] << " in the pedigree";
				throw std::runtime_error(ss.str().c_str()); 
			}
			try
			{
				getFunnel(std::distance(pedigreeLineNames.begin(), i), mother, father, &(funnel.val[0]), nFounders);
			}
			catch(...)
			{
				std::stringstream ss;
				ss << "Attempting to trace pedigree for line " << finalNames(lineCounter) << ": Unable to get funnel";
				throw std::runtime_error(ss.str().c_str());
			}
			allFunnels.push_back(funnel);
		}
		else
		{
			std::vector<long> linesToCheck;
			try
			{
				getAICParentLines(mother, father, pedigreeRow, intercrossingGenerations[lineCounter], linesToCheck);
			}
			catch(...)
			{
				std::stringstream ss;
				ss << "Error while attempting to trace intercrossing for line " << Rcpp::as<std::string>(*i) << " which has " << intercrossingGenerations[lineCounter] << " generations of intercrossing";
				throw std::runtime_error(ss.str().c_str());
			}
			for(std::vector<long>::iterator i = linesToCheck.begin(); i != linesToCheck.end(); i++)
			{
				getFunnel(*i, mother, father, &(funnel.val[0]), nFounders);
				allFunnels.push_back(funnel);
			}

		}
	}
	Rcpp::IntegerMatrix result(allFunnels.size(), nFounders);
	for(std::size_t i = 0; i < allFunnels.size(); i++)
	{
		if(standardise)
		{
			orderFunnel(allFunnels[i].val, nFounders);
		}
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			result(i, founderCounter) = allFunnels[i].val[founderCounter];
		}
	}
	return result;
END_RCPP
}

