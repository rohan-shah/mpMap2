#include "recodeFoundersFinalsHets.h"
#include "testDistortion.h"
SEXP testDistortion(SEXP geneticData_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData = Rcpp::as<Rcpp::S4>(geneticData_sexp);
	Rcpp::IntegerMatrix founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	Rcpp::IntegerMatrix finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	Rcpp::S4 hetData = Rcpp::as<Rcpp::S4>(geneticData.slot("hetData"));
	Rcpp::List hetDataList = Rcpp::as<Rcpp::List>(hetData);
	int nMarkers = hetDataList.size(), nFounders = founders.nrow(), nFinals = finals.nrow();
	if(nMarkers != founders.ncol() || nMarkers != finals.ncol())
	{
		throw std::runtime_error("Number of markers was inconsistent");
	}
	//Testing for distortion is only trivial in the case where there are no hetrozygotes. If there are hetrozygotes then the population may be a mixture due to funnel effects - For certain lines certain hets will be impossible. 
	for(int i = 0; i < nMarkers; i++)
	{
		Rcpp::IntegerMatrix hetDataEntry = Rcpp::as<Rcpp::IntegerMatrix>(hetDataList(i));
		for(int j = 0; j < hetDataEntry.nrow(); j++)
		{
			if(hetDataEntry(j, 0) != hetDataEntry(j, 1)) throw std::runtime_error("Cannot test for segregation distortion if there are hetrozygotes");
		}
	}

	//Re-encoding the data makes things easier
	recodeDataStruct recodeDataArgs;
	recodeDataArgs.hetData = hetData;
	recodeDataArgs.founders = founders;
	recodeDataArgs.finals = finals;
	Rcpp::IntegerMatrix recodedFounders(founders.nrow(), founders.ncol()), recodedFinals(finals.nrow(), finals.ncol());
	Rcpp::List recodedHetData(hetDataList.size());
	recodeDataArgs.recodedFounders = recodedFounders;
	recodeDataArgs.recodedFinals = recodedFinals;
	recodeDataArgs.recodedHetData = recodedHetData;
	try
	{
		recodeFoundersFinalsHets(recodeDataArgs);
	}
	catch(std::invalid_argument& argument)
	{
		throw std::runtime_error("Invalid input, please run validObject on the input mpcross object for more information");
	}

	//Working data
	std::vector<int> alleles, uniqueAlleles, alleleCounts, encodings, table;
	std::vector<double> probabilities;
	
	//Results
	Rcpp::IntegerVector classes(nMarkers);
	Rcpp::NumericVector testStatistics(nMarkers), L1(nMarkers), L2(nMarkers);
	for(int i = 0; i < nMarkers; i++)
	{
		Rcpp::IntegerMatrix hetDataEntry = Rcpp::as<Rcpp::IntegerMatrix>(recodedHetData(i));
		alleles.clear();
		for(int j = 0; j < nFounders; j++)
		{
			alleles.push_back(recodedFounders(j, i));
		}
		uniqueAlleles = alleles;
		std::sort(uniqueAlleles.begin(), uniqueAlleles.end());
		uniqueAlleles.erase(std::unique(uniqueAlleles.begin(), uniqueAlleles.end()), uniqueAlleles.end());
		alleleCounts.resize(uniqueAlleles.size());
		std::fill(alleleCounts.begin(), alleleCounts.end(), 0);
		for(int j = 0; j < nFounders; j++)
		{
			alleleCounts[recodedFounders(j, i)]++;
		}
		probabilities.resize(uniqueAlleles.size());
		for(std::size_t j = 0; j < uniqueAlleles.size(); j++) probabilities[j] = (double)alleleCounts[j] / (double)nFounders;

		table.resize(uniqueAlleles.size());
		std::fill(table.begin(), table.end(), 0);
		int counter = 0;
		for(int j = 0; j < nFinals; j++)
		{
			if(recodedFinals(j, i) != NA_INTEGER)
			{
				table[recodedFinals(j, i)]++;
				counter++;
			}
		}
		testStatistics(i) = L1(i) = L2(i) = 0;
		for(std::size_t j = 0; j < uniqueAlleles.size(); j++)
		{
			double expected = counter*probabilities[j];
			testStatistics(i) += (expected - table[j]) * (expected - table[j]) / expected;
			L1(i) += fabs(probabilities[j] - (double)table[j]/(double)counter);
			double l2Part = probabilities[j] - (double)table[j]/(double)counter;
			L2(i) += l2Part * l2Part;
		}
		classes[i] = uniqueAlleles.size();
		L2(i) = std::sqrt(L2(i));
	}
	return Rcpp::List::create(Rcpp::Named("classes") = classes, Rcpp::Named("testStatistics") = testStatistics, Rcpp::Named("L1") = L1, Rcpp::Named("L2") = L2);
END_RCPP
}
