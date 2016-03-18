#include "multiparentSNP.h"
SEXP multiparentSNPRemoveHets(SEXP object_sexp)
{
BEGIN_RCPP
	Rcpp::S4 object;
	try
	{
		object = Rcpp::as<Rcpp::S4>(object_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input object must be an S4 object");
	}
	if(!object.is("geneticData"))
	{
		throw std::runtime_error("Input object must have class \"geneticData\"");
	}
	Rcpp::S4 copied = Rcpp::clone(object);
	
	Rcpp::IntegerMatrix founders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(copied.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@founders must be an integer matrix");
	}

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(copied.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@finals must be an integer matrix");
	}

	Rcpp::List hetData;
	try
	{
		hetData = Rcpp::as<Rcpp::List>(copied.slot("hetData"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@hetData must be a list");
	}

	Rcpp::Function nMarkersFunc("nMarkers"), sample("sample"), nFoundersFunc("nFounders"), nLinesFunc("nLines");
	int nMarkers = Rcpp::as<int>(nMarkersFunc(copied));
	int nFounders = Rcpp::as<int>(nFoundersFunc(copied));
	int nLines = Rcpp::as<int>(nLinesFunc(copied));

	std::vector<bool> becomesOne(nLines), becomesZero(nLines), isNA(nLines);

	Rcpp::IntegerVector foundersRange = Rcpp::Range(1, nFounders);

	Rcpp::IntegerMatrix newHetData(2, 3);
	newHetData(0, 0) = newHetData(0, 1) = newHetData(0, 2) = 0;
	newHetData(1, 0) = newHetData(1, 1) = newHetData(1, 2) = 1;

	std::vector<int> oneAllelesFounderValues, zeroAllelesFounderValues;
	std::vector<int> oneAllelesFinalValues, zeroAllelesFinalValues;
	for(int i = 0; i < nMarkers; i++)
	{
		int numberOfOnes = Rcpp::as<int>(sample(nFounders-1, 1));
		Rcpp::IntegerVector oneAlleles = Rcpp::as<Rcpp::IntegerVector>(sample(foundersRange, numberOfOnes));
		std::sort(oneAlleles.begin(), oneAlleles.end());
		oneAllelesFounderValues.resize(numberOfOnes);
		zeroAllelesFounderValues.resize(nFounders - numberOfOnes);

		Rcpp::IntegerMatrix relevantHetData = hetData(i);

		//Get out all the old alleles which are going to be converted to ones
		for(int j = 0; j < numberOfOnes; j++) oneAllelesFounderValues[j] = founders(oneAlleles(j)-1, i);
		//Get out all the old alleles which are going to be converted to zeros
		int counter = 0;
		for(int k = 1; k != oneAlleles(0); k++)
		{
			zeroAllelesFounderValues[counter] = k;
			counter++;
		}
		for(int j = 0; j < numberOfOnes-1; j++)
		{
			for(int k = oneAlleles(j) + 1; k != oneAlleles(j+1); k++)
			{
				zeroAllelesFounderValues[counter] = k;
				counter++;
			}
		}
		for(int k = oneAlleles(numberOfOnes-1)+1; k <= nFounders; k++)
		{
			zeroAllelesFounderValues[counter] = k;
			counter++;
		}
		oneAllelesFinalValues.clear();
		zeroAllelesFinalValues.clear();
		//hetValues.clear();
		for(int j = 0; j < relevantHetData.nrow(); j++)
		{
			bool firstZero = std::find(zeroAllelesFounderValues.begin(), zeroAllelesFounderValues.end(), relevantHetData(j, 0)) != zeroAllelesFounderValues.end();
			bool secondZero = std::find(zeroAllelesFounderValues.begin(), zeroAllelesFounderValues.end(), relevantHetData(j, 1)) != zeroAllelesFounderValues.end();
			if(firstZero && secondZero) zeroAllelesFinalValues.push_back(relevantHetData(j, 2));
			//else if(firstZero ^ secondZero) hetValues.push_back(relevantHetData(j, 2));
			else if(!(firstZero ^ secondZero)) oneAllelesFinalValues.push_back(relevantHetData(j, 2));
		}

		std::fill(becomesOne.begin(), becomesOne.end(), false);
		std::fill(becomesZero.begin(), becomesZero.end(), false);
		std::fill(isNA.begin(), isNA.end(), false);
		//Work out which lines have an NA value
		for(int k = 0; k < nLines; k++)
		{
			if(finals(k, i) == NA_INTEGER) isNA[k] = true;
		}
		//Work out which of the finals become 1
		for(int j = 0; j < (int)oneAllelesFinalValues.size(); j++)
		{
			for(int k = 0; k < nLines; k++)
			{
				if(finals(k, i) == oneAllelesFinalValues[j])
				{
					becomesOne[k] = true;
				}
			}
		}
		//..and which become 0
		for(int j = 0; j < (int)zeroAllelesFinalValues.size(); j++)
		{
			for(int k = 0; k < nLines; k++)
			{
				if(finals(k, i) == zeroAllelesFinalValues[j])
				{
					becomesZero[k] = true;
				}
			}
		}
		//The remainder will become NA
		for(int j = 0; j < nLines; j++) 
		{
			if(isNA[j]) finals(j, i) = NA_INTEGER;
			else if(becomesOne[j]) finals(j, i) = 1;
			else if(becomesZero[j]) finals(j, i) = 0;
			//This is the het case
			else finals(j, i) = NA_INTEGER;
		}
		//Fill the founders with zeros....
		for(int j = 0; j < nFounders; j++) founders(j, i) = 0;
		//...and put in the ones
		for(int j = 0; j < oneAlleles.size(); j++) founders(oneAlleles[j]-1, i) = 1;

		hetData(i) = newHetData;
	}
	return copied;
END_RCPP
}
SEXP multiparentSNPKeepHets(SEXP object_sexp)
{
BEGIN_RCPP
	Rcpp::S4 object;
	try
	{
		object = Rcpp::as<Rcpp::S4>(object_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input object must be an S4 object");
	}
	if(!object.is("geneticData"))
	{
		throw std::runtime_error("Input object must have class \"geneticData\"");
	}
	Rcpp::S4 copied = Rcpp::clone(object);
	
	Rcpp::IntegerMatrix founders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(copied.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@founders must be an integer matrix");
	}

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(copied.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@finals must be an integer matrix");
	}

	Rcpp::List hetData;
	try
	{
		hetData = Rcpp::as<Rcpp::List>(copied.slot("hetData"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot object@hetData must be a list");
	}

	Rcpp::Function nMarkersFunc("nMarkers"), sample("sample"), nFoundersFunc("nFounders"), nLinesFunc("nLines");
	int nMarkers = Rcpp::as<int>(nMarkersFunc(copied));
	int nFounders = Rcpp::as<int>(nFoundersFunc(copied));
	int nLines = Rcpp::as<int>(nLinesFunc(copied));

	std::vector<bool> becomesOne(nLines), becomesZero(nLines), isNA(nLines);

	Rcpp::IntegerVector foundersRange = Rcpp::Range(1, nFounders);

	Rcpp::IntegerMatrix newHetData(4, 3);
	newHetData(0, 0) = newHetData(0, 1) = newHetData(0, 2) = 0;
	newHetData(1, 0) = newHetData(1, 1) = newHetData(1, 2) = 1;
	newHetData(2, 0) = 0; newHetData(2, 1) = 1; newHetData(2, 2) = 2;
	newHetData(3, 0) = 1; newHetData(3, 1) = 0; newHetData(3, 2) = 2;

	std::vector<int> oneAllelesFounderValues, zeroAllelesFounderValues;
	std::vector<int> oneAllelesFinalValues, zeroAllelesFinalValues;
	//std::vector<int> hetValues;
	for(int i = 0; i < nMarkers; i++)
	{
		int numberOfOnes = Rcpp::as<int>(sample(nFounders-1, 1));
		Rcpp::IntegerVector oneAlleles = Rcpp::as<Rcpp::IntegerVector>(sample(foundersRange, numberOfOnes));
		std::sort(oneAlleles.begin(), oneAlleles.end());
		oneAllelesFounderValues.resize(numberOfOnes);
		zeroAllelesFounderValues.resize(nFounders - numberOfOnes);

		Rcpp::IntegerMatrix relevantHetData = hetData(i);

		//Get out all the old alleles which are going to be converted to ones
		for(int j = 0; j < numberOfOnes; j++) oneAllelesFounderValues[j] = founders(oneAlleles(j)-1, i);
		//Get out all the old alleles which are going to be converted to zeros
		int counter = 0;
		for(int k = 1; k != oneAlleles(0); k++)
		{
			zeroAllelesFounderValues[counter] = k;
			counter++;
		}
		for(int j = 0; j < numberOfOnes-1; j++)
		{
			for(int k = oneAlleles(j) + 1; k != oneAlleles(j+1); k++)
			{
				zeroAllelesFounderValues[counter] = k;
				counter++;
			}
		}
		for(int k = oneAlleles(numberOfOnes-1)+1; k <= nFounders; k++)
		{
			zeroAllelesFounderValues[counter] = k;
			counter++;
		}
		oneAllelesFinalValues.clear();
		zeroAllelesFinalValues.clear();
		//hetValues.clear();
		for(int j = 0; j < relevantHetData.nrow(); j++)
		{
			bool firstZero = std::find(zeroAllelesFounderValues.begin(), zeroAllelesFounderValues.end(), relevantHetData(j, 0)) != zeroAllelesFounderValues.end();
			bool secondZero = std::find(zeroAllelesFounderValues.begin(), zeroAllelesFounderValues.end(), relevantHetData(j, 1)) != zeroAllelesFounderValues.end();
			if(firstZero && secondZero) zeroAllelesFinalValues.push_back(relevantHetData(j, 2));
			//else if(firstZero ^ secondZero) hetValues.push_back(relevantHetData(j, 2));
			else if(!(firstZero ^ secondZero)) oneAllelesFinalValues.push_back(relevantHetData(j, 2));
		}
		std::fill(becomesOne.begin(), becomesOne.end(), false);
		std::fill(becomesZero.begin(), becomesZero.end(), false);
		std::fill(isNA.begin(), isNA.end(), false);
		//Work out which lines have an NA value
		for(int k = 0; k < nLines; k++)
		{
			if(finals(k, i) == NA_INTEGER) isNA[k] = true;
		}
		//Work out which of the finals become 1
		for(int j = 0; j < (int)oneAllelesFinalValues.size(); j++)
		{
			for(int k = 0; k < nLines; k++)
			{
				if(finals(k, i) == oneAllelesFinalValues[j])
				{
					becomesOne[k] = true;
				}
			}
		}
		//..and which become 0
		for(int j = 0; j < (int)zeroAllelesFinalValues.size(); j++)
		{
			for(int k = 0; k < nLines; k++)
			{
				if(finals(k, i) == zeroAllelesFinalValues[j])
				{
					becomesZero[k] = true;
				}
			}
		}
		//The remainder will become NA
		for(int j = 0; j < nLines; j++) 
		{
			if(isNA[j]) finals(j, i) = NA_INTEGER;
			else if(becomesOne[j]) finals(j, i) = 1;
			else if(becomesZero[j]) finals(j, i) = 0;
			//This is the het case
			else finals(j, i) = 2;
		}
		//Fill the founders with zeros....
		for(int j = 0; j < nFounders; j++) founders(j, i) = 0;
		//...and put in the ones
		for(int j = 0; j < oneAlleles.size(); j++) founders(oneAlleles[j]-1, i) = 1;

		hetData(i) = newHetData;
	}
	return copied;
END_RCPP
}


