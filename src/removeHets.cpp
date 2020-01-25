#include "removeHets.h"
#include "recodeFoundersFinalsHets.h"
SEXP removeHets(SEXP founders_sexp, SEXP finals_sexp, SEXP hetData_sexp)
{
BEGIN_RCPP
	Rcpp::IntegerMatrix founders = Rcpp::as<Rcpp::IntegerMatrix>(founders_sexp);
	Rcpp::IntegerMatrix finals = Rcpp::as<Rcpp::IntegerMatrix>(finals_sexp);
	Rcpp::S4 hetData = Rcpp::as<Rcpp::S4>(hetData_sexp);
	Rcpp::List hetDataList = Rcpp::as<Rcpp::List>(hetData_sexp);
	int nMarkers = founders.ncol(), nFounders = founders.nrow(), nFinals = finals.nrow();

	recodeDataStruct args;
	Rcpp::List recodedHetData(nMarkers);
	Rcpp::IntegerMatrix recodedFounders(founders.nrow(), founders.ncol()), recodedFinals(finals.nrow(), finals.ncol());
	args.founders = founders;
	args.finals = finals;
	args.hetData = hetData;
	args.recodedHetData = recodedHetData;
	args.recodedFounders = recodedFounders;
	args.recodedFinals = recodedFinals;
	try
	{
		recodeFoundersFinalsHets(args);
	}
	catch(std::invalid_argument& argument)
	{
		throw std::runtime_error("Invalid input, please run validObject on the input mpcross object for more information");
	}

	for(int i = 0; i < nMarkers; i++)
	{
		int maxAllele = 0;
		for(int j = 0; j < nFounders; j++) maxAllele = std::max(maxAllele, recodedFounders(j, i));
		//Overwrite the recoded finals data
		for(int j = 0; j < nFinals; j++)
		{
			if(recodedFinals(j, i) > maxAllele) recodedFinals(j, i) = NA_INTEGER;
			else recodedFinals(j, i) = finals(j, i);
		}
		//overwrite the recoded het data
		Rcpp::IntegerMatrix newHetEntry(maxAllele+1, 3);
		Rcpp::IntegerMatrix oldHetEntry = Rcpp::as<Rcpp::IntegerMatrix>(hetDataList(i));
		int counter = 0;
		for(int j = 0; j < oldHetEntry.nrow(); j++)
		{
			if(oldHetEntry(j, 0) == oldHetEntry(j, 1))
			{
				newHetEntry(counter, 0) = newHetEntry(counter, 1) = newHetEntry(counter, 2) = oldHetEntry(j, 1);
				counter++;
			}
		}
		recodedHetData(i) = newHetEntry;
	}
	Rcpp::rownames(recodedFinals) = Rcpp::rownames(finals);
	Rcpp::colnames(recodedFinals) = Rcpp::colnames(finals);
	return Rcpp::List::create(Rcpp::Named("finals") = recodedFinals, Rcpp::Named("hetData") = recodedHetData);
END_RCPP
}
