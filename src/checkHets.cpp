#include "checkHets.h"
SEXP checkHets(SEXP hets)
{
	std::vector<std::string> errors;
	Rcpp::List hetObject = hets;
	Rcpp::CharacterVector hetObjectNames = hetObject.names();
	std::size_t index = 0;
	for(Rcpp::List::iterator i = hetObject.begin(); i != hetObject.end(); i++,index++)
	{
		if(errors.size() > 10) break;
		Rcpp::RObject currentHetObject = *i;
		if(currentHetObject.sexp_type() != INTSXP || !currentHetObject.hasAttribute("dim"))
		{
			errors.push_back("Entry for marker " + Rcpp::as<std::string>(hetObjectNames[index]) + " was not an integer matrix");
			continue;
		}
		Rcpp::IntegerVector currentHetObjectDim = Rcpp::as<Rcpp::IntegerVector>(currentHetObject.attr("dim"));
		if(currentHetObjectDim.size() != 2 || currentHetObjectDim[1] != 3)
		{
			errors.push_back("Entry for marker " + Rcpp::as<std::string>(hetObjectNames[index]) + " must be a matrix with three columns");
			continue;
		}
		Rcpp::IntegerMatrix currentHetObjectMat = Rcpp::as<Rcpp::IntegerMatrix>(currentHetObject);
		int currentHetObjectRows = currentHetObjectDim(0);
		//Now check symmetry - That is, if we have an encoding for haplotype (a, b), then we must have the same encoding for haplotype (b, a)
		for(int rowCounter = 0; rowCounter < currentHetObjectRows; rowCounter++)
		{
			int otherRow = 0;
			for(; otherRow < currentHetObjectRows; otherRow++)
			{
				if(currentHetObjectMat(otherRow, 0) == currentHetObjectMat(rowCounter, 1) && currentHetObjectMat(otherRow, 1) == currentHetObjectMat(rowCounter, 0) && currentHetObjectMat(otherRow, 2) == currentHetObjectMat(rowCounter, 2))
				{
					break;
				}
			}
			//If we didn't find the symmetric encoding, that's an error. 
			if(otherRow == currentHetObjectRows) 
			{
				errors.push_back("Entry for marker " + Rcpp::as<std::string>(hetObjectNames[index]) + ": If haplotype (a, b) has an encoding then haplotype (b, a) must have the same encoding");
				break;
			}

		}
	}
	if(errors.size() > 0)
	{
		Rcpp::CharacterVector retVal = Rcpp::wrap(errors);
		return retVal;
	}
	return Rcpp::wrap(true);
}