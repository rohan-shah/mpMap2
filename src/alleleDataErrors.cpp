#include "alleleDataErrors.h"
#include <vector>
#include <string>
RcppExport SEXP alleleDataErrors(SEXP Robject, SEXP Rlimit)
{
	BEGIN_RCPP
		Rcpp::S4 mpcross = Robject;
		Rcpp::IntegerMatrix founders = mpcross.slot("founders"), finals = mpcross.slot("finals");
		Rcpp::List hetData = mpcross.slot("hetData");
		Rcpp::CharacterVector markerNames = hetData.attr("names");
		int limit = Rcpp::as<int>(Rlimit);

		Rcpp::List codingErrors = listCodingErrors(founders, finals, hetData);
		Rcpp::IntegerMatrix founderErrors = codingErrors["founders"];
		Rcpp::IntegerMatrix finalErrors = codingErrors["finals"];

		std::vector<std::string> codingErrorsAsStrings;
		int nFounderErrors = founderErrors.nrow(), nFinalErrors = finalErrors.nrow();
		for(int i = 0; i < nFounderErrors; i++)
		{
			if((int)codingErrorsAsStrings.size() >= limit) 
			{
				codingErrorsAsStrings.push_back("Omitting details of further coding errors");
				break;
			}
			int markerIndex = founderErrors(i, 0), hetDataRow = founderErrors(i, 1), hetDataColumn = founderErrors(i, 2);
			Rcpp::IntegerMatrix hetDataEntry = hetData(markerIndex);
			std::stringstream ss;
			ss << "Coding error for marker " << markerNames[markerIndex] << ": Value " << hetDataEntry(hetDataRow, hetDataColumn) << " of @hetData[[" << markerIndex+1 << "]][" << hetDataRow+1 << ", " << hetDataColumn+1 << "] not present in @founders[,"<< markerIndex+1<<"]";
			codingErrorsAsStrings.push_back(ss.str());
		}
		if((int)codingErrorsAsStrings.size() < limit)
		{
			for(int i = 0; i < nFinalErrors; i++)
			{
				if((int)codingErrorsAsStrings.size() >= limit) 
				{
					codingErrorsAsStrings.push_back("Omitting details of further coding errors");
					break;
				}
				Rcpp::IntegerMatrix hetDataEntry = hetData(finalErrors(i, 1));
				std::stringstream ss;
				ss << "Coding error for marker " << markerNames[finalErrors(i, 1)] << ": Value " << finals(finalErrors(i, 0), finalErrors(i, 1)) << " from finals slot not present in hetData slot";
				codingErrorsAsStrings.push_back(ss.str());
			}
		}
		return Rcpp::wrap(codingErrorsAsStrings);
	END_RCPP
}
RcppExport SEXP listCodingErrors(SEXP _founders, SEXP _finals, SEXP _hetData)
{
	BEGIN_RCPP
		Rcpp::IntegerMatrix founders = _founders, finals = _finals;
		Rcpp::List hetData = _hetData;

		std::vector<int> founderErrorMarker, founderErrorRow, founderErrorColumn, finalErrorRow, finalErrorMarker;
		int nMarkers = hetData.length();
		int nFinals = finals.nrow(), nFounders = founders.nrow();
		//Possible values for finals
		std::vector<int> validFinalValues, validFounderValues;
		for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
		{
			validFounderValues.clear();
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				validFounderValues.push_back(founders(founderCounter, markerCounter));
			}
			std::sort(validFounderValues.begin(), validFounderValues.end());
			validFounderValues.erase(std::unique(validFounderValues.begin(), validFounderValues.end()), validFounderValues.end());
			Rcpp::IntegerMatrix currentHetData = hetData(markerCounter);

			validFinalValues.clear();
			for(int hetRowCounter = 0; hetRowCounter < currentHetData.nrow(); hetRowCounter++)
			{
				if(std::find(validFounderValues.begin(), validFounderValues.end(), currentHetData(hetRowCounter, 0)) == validFounderValues.end())
				{
					founderErrorMarker.push_back(markerCounter);
					founderErrorRow.push_back(hetRowCounter);
					founderErrorColumn.push_back(0);
				}
				if(std::find(validFounderValues.begin(), validFounderValues.end(), currentHetData(hetRowCounter, 1)) == validFounderValues.end())
				{
					founderErrorMarker.push_back(markerCounter);
					founderErrorRow.push_back(hetRowCounter);
					founderErrorColumn.push_back(1);
				}
				validFinalValues.push_back(currentHetData(hetRowCounter, 2));
			}
			for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
			{
				if(std::find(validFinalValues.begin(), validFinalValues.end(), finals(finalCounter, markerCounter)) == validFinalValues.end())
				{
					finalErrorRow.push_back(finalCounter);
					finalErrorMarker.push_back(markerCounter);
				}
			}
		}
		Rcpp::IntegerMatrix founderErrors(founderErrorRow.size(), 3), finalErrors(finalErrorRow.size(), 2);
		for(std::size_t founderErrorCounter = 0; founderErrorCounter < founderErrorRow.size(); founderErrorCounter++)
		{
			founderErrors(founderErrorCounter, 0) = founderErrorMarker[founderErrorCounter];
			founderErrors(founderErrorCounter, 1) = founderErrorRow[founderErrorCounter];
			founderErrors(founderErrorCounter, 2) = founderErrorColumn[founderErrorCounter];
		}
		founderErrors.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("Marker", "Row", "Column"));

		for(std::size_t i = 0; i < finalErrorRow.size(); i++)
		{
			finalErrors(i, 0) = finalErrorRow[i];
			finalErrors(i, 1) = finalErrorMarker[i];
		}
		finalErrors.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("Row", "Column"));
		return Rcpp::List::create(Rcpp::Named("founders") = founderErrors, Rcpp::Named("finals") = finalErrors);
	END_RCPP
}