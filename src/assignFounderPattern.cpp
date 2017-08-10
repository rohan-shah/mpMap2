#include "assignFounderPattern.h"
RcppExport SEXP assignFounderPattern(SEXP geneticData_sexp, SEXP newFounders_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData;
	try
	{
		geneticData = Rcpp::as<Rcpp::S4>(geneticData_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData must be an S4 object");
	}

	Rcpp::IntegerMatrix founders, finals, newFounders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData must have a slot named founders, containing an integer matrix");
	}

	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData must have a slot named finals, containing an integer matrix");
	}

	try
	{
		newFounders = Rcpp::as<Rcpp::IntegerMatrix>(newFounders_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input newFounders must be an integer matrix");
	}

	Rcpp::List hetData;
	try
	{
		hetData = Rcpp::as<Rcpp::List>(geneticData.slot("hetData"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@hetData must be a list");
	}
	int nFounders = founders.nrow();
	int maxAlleles = nFounders*(nFounders+1)/2;

	int newFoundersRows = newFounders.nrow(), newFoundersCols = newFounders.ncol(), foundersCols = newFounders.ncol();
	if(newFoundersRows != nFounders || newFoundersCols != foundersCols)
	{
		throw std::runtime_error("Input newFounders had the wrong dimensions");
	}

	for(int marker = 0; marker < newFoundersCols; marker++)
	{
		for(int founder = 0; founder < nFounders; founder++)
		{
			if(founders(founder, marker) != founder+1) throw std::runtime_error("Can only apply assignFounderPattern to an object with fully informative founders");
		}
	}

	Rcpp::CharacterVector markerNames = Rcpp::colnames(finals);

	Rcpp::IntegerMatrix newFinals(finals.nrow(), finals.ncol());
	Rcpp::List newHetDataList(finals.ncol());
	newHetDataList.names() = markerNames;

	std::vector<int> translationVector(maxAlleles);
	std::vector<int> allAlleles, newHomozygotes;
	for(int marker = 0; marker < foundersCols; marker++)
	{
		Rcpp::IntegerMatrix currentHetData = Rcpp::as<Rcpp::IntegerMatrix>(hetData(marker));
		std::fill(translationVector.begin(), translationVector.end(), -1);
		for(int founder = 0; founder < nFounders; founder++) translationVector[founder] = newFounders(founder, marker);

		//Get out a list of all the alleles (originally).
		allAlleles.resize(currentHetData.nrow());
		for(int i = 0; i < currentHetData.nrow(); i++) allAlleles[i] = currentHetData(i, 2);
		std::sort(allAlleles.begin(), allAlleles.end());
		allAlleles.erase(std::unique(allAlleles.begin(), allAlleles.end()), allAlleles.end());

		//Check that they're all in the expected range.
		for(std::size_t i = 0; i < allAlleles.size(); i++)
		{
			if(allAlleles[i] < 1 || allAlleles[i] > maxAlleles)
			{
				throw std::runtime_error("An encoding for an allele was out of range. All alleles should be between 0 and (nFounders*(nFounders+1)/2 + 1)");
			}
		}

		newHomozygotes.resize(nFounders);
		for(int founder = 0; founder < nFounders; founder++) newHomozygotes[founder] = newFounders(founder, marker);
		std::sort(newHomozygotes.begin(), newHomozygotes.end());
		newHomozygotes.erase(std::unique(newHomozygotes.begin(), newHomozygotes.end()), newHomozygotes.end());
		int newHet = *std::max_element(newHomozygotes.begin(), newHomozygotes.end())+1;

		Rcpp::IntegerMatrix currentNewHetData((allAlleles.size() - nFounders)*2 + nFounders, 3);
		for(int i = 0; i < currentHetData.nrow(); i++)
		{
			currentNewHetData(i, 0) = translationVector[currentHetData(i, 0)-1];
			currentNewHetData(i, 1) = translationVector[currentHetData(i, 1)-1];
		}

		for(int i = 0; i < currentHetData.nrow(); i++)
		{
			if(translationVector[currentHetData(i, 2) - 1] == -1)
			{
				int firstRow = 0;
				for(; firstRow <= i; firstRow++)
				{
					if(currentNewHetData(firstRow, 0) == currentNewHetData(i, 0) && currentNewHetData(firstRow, 1) == currentNewHetData(i, 1))
					{
						break;
					}
				}
				if(translationVector[currentHetData(firstRow, 2) - 1] == -1)
				{
					translationVector[currentHetData(i, 2) - 1] = newHet;
					newHet++;
				}
				else
				{
					translationVector[currentHetData(i, 2) - 1] = translationVector[currentHetData(firstRow, 2) - 1];
				}
			}
		}
		for(int i = 0; i < currentNewHetData.nrow(); i++) currentNewHetData(i, 2) = translationVector[currentHetData(i, 2) - 1];
		int destinationRow = 0;
		for(int i = 0; i < currentNewHetData.nrow(); i++)
		{
			bool keep = true;
			for(int j = 0; j < i; j++)
			{
				if(currentNewHetData(j, 0) == currentNewHetData(i, 0) && currentNewHetData(j, 1) == currentNewHetData(i, 1))
				{
					keep = false;
					break;
				}
			}
			if(keep)
			{
				currentNewHetData(destinationRow, 0) = currentNewHetData(i, 0);
				currentNewHetData(destinationRow, 1) = currentNewHetData(i, 1);
				currentNewHetData(destinationRow, 2) = currentNewHetData(i, 2);
				destinationRow++;
			}
		}
		newHetDataList(marker) = currentNewHetData(Rcpp::Range(0, destinationRow-1), Rcpp::_);

		for(int i = 0; i < finals.nrow(); i++) newFinals(i, marker) = translationVector[finals(i, marker) - 1];
	}
	Rcpp::S4 newHetData("hetData");
	newHetData.slot(".Data") = newHetDataList;

	Rcpp::S4 newDataObject("geneticData");

	Rcpp::rownames(newFinals) = Rcpp::rownames(finals);
	Rcpp::colnames(newFinals) = Rcpp::colnames(finals);
	newDataObject.slot("finals") = newFinals;

	newDataObject.slot("founders") = newFounders;
	newDataObject.slot("pedigree") = geneticData.slot("pedigree");
	newDataObject.slot("hetData") = newHetData;
	return newDataObject;
END_RCPP
}
