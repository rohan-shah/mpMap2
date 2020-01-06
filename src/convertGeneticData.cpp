#include "convertGeneticData.h"
void convertGeneticData(Rcpp::S4 obj)
{
	Rcpp::RObject finals = obj.slot("finals");
	Rcpp::Function convert("storage.mode<-");
	if(finals.sexp_type() == REALSXP)
	{
		finals = obj.slot("finals") = convert(finals, "integer");
	}
	Rcpp::RObject founders = obj.slot("founders");
	if(founders.sexp_type() == REALSXP)
	{
		founders = obj.slot("founders") = convert(founders, "integer");
	}
	Rcpp::List hetData = obj.slot("hetData");
	for(int i = 0; i < hetData.size(); i++)
	{
		if(Rcpp::as<Rcpp::RObject>(hetData(i)).sexp_type() == REALSXP)
		{
			hetData(i) = convert(Rcpp::as<Rcpp::NumericMatrix>(hetData(i)), "integer");
		}
	}
}

