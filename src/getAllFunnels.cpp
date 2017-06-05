#include "getAllFunnels.h"
#include "intercrossingAndSelfingGenerations.h"
#include "estimateRF.h"
#include "getFunnel.h"
#include "orderFunnel.h"
SEXP getAllFunnels(SEXP geneticData_sexp, SEXP standardise_sexp)
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

	int nFounders = founders.nrow(), nFinals = finals.nrow();
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(pedigree, finals, nFounders, intercrossingGenerations, selfingGenerations);

	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("mother"));
	Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(pedigree.slot("father"));

	funnelType funnel;
	Rcpp::IntegerMatrix result(nFinals, nFounders);

	for(int lineCounter = 0; lineCounter < finalNames.size(); lineCounter++)
	{
		if(intercrossingGenerations[lineCounter] == 0)
		{
			Rcpp::CharacterVector::iterator i = std::find(pedigreeLineNames.begin(), pedigreeLineNames.end(), finalNames[lineCounter]);
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
			if(standardise)
			{
				orderFunnel(funnel.val, nFounders);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++) result(lineCounter, founderCounter) = funnel.val[founderCounter];
		}
		else
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++) result(lineCounter, founderCounter) = NA_INTEGER;
		}
	}
	return result;
END_RCPP
}

