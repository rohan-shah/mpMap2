#include "alleleDataErrors.h"
#include <vector>
#include <string>
#include "convertGeneticData.h"
RcppExport SEXP alleleDataErrors(SEXP Robject, SEXP Rlimit)
{
	BEGIN_RCPP
		Rcpp::S4 geneticData = Robject;
		convertGeneticData(geneticData);

		Rcpp::IntegerMatrix founders = geneticData.slot("founders"), finals = geneticData.slot("finals");
		Rcpp::List hetData = geneticData.slot("hetData");
		int limit = Rcpp::as<int>(Rlimit);

		Rcpp::List codingErrors = listCodingErrors(founders, finals, hetData);
		std::vector<std::string> codingErrorsAsStrings;
		codingErrorsToStrings(codingErrors, codingErrorsAsStrings, finals, hetData, limit);
		return Rcpp::wrap(codingErrorsAsStrings);
	END_RCPP
}
RcppExport void codingErrorsToStrings(Rcpp::List codingErrors, std::vector<std::string>& codingErrorsAsStrings, Rcpp::IntegerMatrix finals, Rcpp::List hetData, int limit)
{
	BEGIN_RCPP
		Rcpp::CharacterVector markerNames = hetData.attr("names");

		Rcpp::IntegerMatrix hetDataErrors = codingErrors["hetData"];
		Rcpp::IntegerMatrix finalErrors = codingErrors["finals"];
		Rcpp::IntegerVector nullErrors = codingErrors["null"];

		R_xlen_t nFounderErrors = hetDataErrors.nrow(), nFinalErrors = finalErrors.nrow(), nNullErrors = nullErrors.length();
		for(int i = 0; i < nFounderErrors; i++)
		{
			if((int)codingErrorsAsStrings.size() >= limit) 
			{
				codingErrorsAsStrings.push_back("Omitting details of further coding errors");
				return;
			}
			int markerIndex = hetDataErrors(i, 0), hetDataRow = hetDataErrors(i, 1), hetDataColumn = hetDataErrors(i, 2);
			Rcpp::IntegerMatrix hetDataEntry = hetData(markerIndex);
			std::stringstream ss;
			ss << "Coding error for marker " << markerNames[markerIndex] << ": Value " << hetDataEntry(hetDataRow, hetDataColumn) << " of @hetData[[" << markerIndex+1 << "]][" << hetDataRow+1 << ", " << hetDataColumn+1 << "] not present in @founders[,"<< markerIndex+1<<"]";
			codingErrorsAsStrings.push_back(ss.str());
		}
		for(int i = 0; i < nFinalErrors; i++)
		{
			if((int)codingErrorsAsStrings.size() >= limit) 
			{
				codingErrorsAsStrings.push_back("Omitting details of further coding errors");
				return;
			}
			Rcpp::IntegerMatrix hetDataEntry = hetData(finalErrors(i, 1));
			std::stringstream ss;
			ss << "Coding error for marker " << markerNames[finalErrors(i, 1)] << ": Value " << finals(finalErrors(i, 0), finalErrors(i, 1)) << " from finals slot not present in hetData slot";
			codingErrorsAsStrings.push_back(ss.str());
		}
		for(int i = 0; i < nNullErrors; i++)
		{
			if((int)codingErrorsAsStrings.size() >= limit) 
			{
				codingErrorsAsStrings.push_back("Omitting details of further coding errors");
				return;
			}
			std::stringstream ss;
			ss << "Coding error for marker " << markerNames[nullErrors(i)] << ": If a founder allele is coded as NA, all founder alleles must be coded as NA, and slot hetData[[" << markerNames[nullErrors(i)] << "]] must be a 0 x 3 matrix";
			codingErrorsAsStrings.push_back(ss.str());
		}
	VOID_END_RCPP
}
RcppExport SEXP listCodingErrors(SEXP _founders, SEXP _finals, SEXP _hetData)
{
	BEGIN_RCPP
		Rcpp::IntegerMatrix founders = _founders, finals = _finals;
		Rcpp::List hetData = _hetData;

		std::vector<int> hetDataErrorMarker, hetDataErrorRow, hetDataErrorColumn;
		std::vector<int> finalErrorRow, finalErrorMarker;
		std::vector<int> nullErrorMarkers;

		R_xlen_t nMarkers = hetData.length();
		int nFinals = finals.nrow(), nFounders = founders.nrow();
		//Possible values for finals
		std::vector<int> validFinalValues, validFounderValues;
		for(R_xlen_t markerCounter = 0; markerCounter < nMarkers; markerCounter++)
		{
			bool hasNullFounder = false;
			validFounderValues.clear();
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				validFounderValues.push_back(founders(founderCounter, markerCounter));
				hasNullFounder |= founders(founderCounter, markerCounter) == NA_INTEGER;
			}
			std::sort(validFounderValues.begin(), validFounderValues.end());
			validFounderValues.erase(std::unique(validFounderValues.begin(), validFounderValues.end()), validFounderValues.end());
			Rcpp::IntegerMatrix currentHetData = hetData(markerCounter);

			validFinalValues.clear();
			if(hasNullFounder)
			{
				if(currentHetData.nrow() != 0)
				{
					nullErrorMarkers.push_back((int)markerCounter);
				}
				else
				{
					for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
					{
						if(finals(finalCounter, markerCounter) != NA_INTEGER)
						{
							nullErrorMarkers.push_back((int)markerCounter);
							break;
						}
					}
				}
			}
			else
			{
				for(int hetRowCounter = 0; hetRowCounter < currentHetData.nrow(); hetRowCounter++)
				{
					if(std::find(validFounderValues.begin(), validFounderValues.end(), currentHetData(hetRowCounter, 0)) == validFounderValues.end())
					{
						hetDataErrorMarker.push_back((int)markerCounter);
						hetDataErrorRow.push_back(hetRowCounter);
						hetDataErrorColumn.push_back(0);
					}
					if(std::find(validFounderValues.begin(), validFounderValues.end(), currentHetData(hetRowCounter, 1)) == validFounderValues.end())
					{
						hetDataErrorMarker.push_back((int)markerCounter);
						hetDataErrorRow.push_back(hetRowCounter);
						hetDataErrorColumn.push_back(1);
					}
					validFinalValues.push_back(currentHetData(hetRowCounter, 2));
				}
				for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
				{
					if(finals(finalCounter, markerCounter) != NA_INTEGER && std::find(validFinalValues.begin(), validFinalValues.end(), finals(finalCounter, markerCounter)) == validFinalValues.end())
					{
						finalErrorRow.push_back(finalCounter);
						finalErrorMarker.push_back((int)markerCounter);
					}
				}
			}
		}
		Rcpp::IntegerMatrix hetDataErrors((int)hetDataErrorRow.size(), 3), finalErrors((int)finalErrorRow.size(), 2);
		Rcpp::IntegerVector nullErrors = Rcpp::wrap(nullErrorMarkers);
		for(int hetDataErrorCounter = 0; hetDataErrorCounter < (int)hetDataErrorRow.size(); hetDataErrorCounter++)
		{
			hetDataErrors(hetDataErrorCounter, 0) = hetDataErrorMarker[hetDataErrorCounter];
			hetDataErrors(hetDataErrorCounter, 1) = hetDataErrorRow[hetDataErrorCounter];
			hetDataErrors(hetDataErrorCounter, 2) = hetDataErrorColumn[hetDataErrorCounter];
		}
		hetDataErrors.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("Marker", "Row", "Column"));

		for(int i = 0; i < (int)finalErrorRow.size(); i++)
		{
			finalErrors(i, 0) = finalErrorRow[i];
			finalErrors(i, 1) = finalErrorMarker[i];
		}
		finalErrors.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("Row", "Column"));
		return Rcpp::List::create(Rcpp::Named("hetData") = hetDataErrors, Rcpp::Named("finals") = finalErrors, Rcpp::Named("null") = nullErrors);
	END_RCPP
}
