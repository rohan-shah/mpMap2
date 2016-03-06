#include "imputeFounders.h"
#include "intercrossingAndSelfingGenerations.h"
#include "recodeFoundersFinalsHets.h"
#include "matrices.hpp"
#include "probabilities.hpp"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#include "funnelsToUniqueValues.h"
#include "estimateRFCheckFunnels.h"
#include "markerPatternsToUniqueValues.h"
#include "intercrossingHaplotypeToMarker.hpp"
#include "funnelHaplotypeToMarker.hpp"
#include "viterbi.hpp"
template<int nFounders, bool infiniteSelfing> void imputedFoundersInternal2(Rcpp::IntegerMatrix founders, Rcpp::IntegerMatrix finals, Rcpp::S4 pedigree, Rcpp::List hetData, Rcpp::List map, Rcpp::IntegerMatrix results)
{
	//Work out maximum number of markers per chromosome
	int maxChromosomeMarkers = 0;
	for(int i = 0; i < map.size(); i++)
	{
		Rcpp::NumericVector chromosome = map(i);
		maxChromosomeMarkers = std::max((int)chromosome.size(), maxChromosomeMarkers);
	}

	typedef typename expandedProbabilities<nFounders, infiniteSelfing>::type expandedProbabilitiesType;
	//expandedProbabilitiesType haplotypeProbabilities;

	Rcpp::Function diff("diff"), haldaneToRf("haldaneToRf");

	//Get out generations of selfing and intercrossing
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(pedigree, finals, nFounders, intercrossingGenerations, selfingGenerations);

	int maxSelfing = *std::max_element(selfingGenerations.begin(), selfingGenerations.end());
	int minSelfing = *std::min_element(selfingGenerations.begin(), selfingGenerations.end());
	int maxAIGenerations = *std::max_element(intercrossingGenerations.begin(), intercrossingGenerations.end());
	int minAIGenerations = *std::min_element(intercrossingGenerations.begin(), intercrossingGenerations.end());
	int nMarkers = founders.ncol();
	int nFinals = finals.nrow();

	//re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles
	//We do this to make it easier to identify markers with identical segregation patterns. recodedFounders = column major matrix
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	Rcpp::List recodedHetData(nMarkers);
	recodedHetData.attr("names") = hetData.attr("names");
	recodedFinals.attr("dimnames") = finals.attr("dimnames");

	recodeDataStruct recoded;
	recoded.recodedFounders = recodedFounders;
	recoded.recodedFinals = recodedFinals;
	recoded.founders = founders;
	recoded.finals = finals;
	recoded.hetData = hetData;
	recoded.recodedHetData = recodedHetData;
	recodeFoundersFinalsHets(recoded);

	//Get out the number of unique funnels. This is only needed because in the case of one funnel we assume a single funnel design and in the case of multiple funnels we assume a random funnels design
	std::vector<std::string> errors, warnings;
	std::vector<funnelType> allFunnels, lineFunnels;
	{
		estimateRFCheckFunnels(recodedFinals, recodedFounders, recodedHetData, pedigree, intercrossingGenerations, warnings, errors, allFunnels, lineFunnels);
		if(errors.size() > 0)
		{
			throw std::runtime_error(errors[0].c_str());
		}
		//Don't bother outputting warnings here
	}
	//map containing encodings of the funnels involved in the experiment (as key), and an associated unique index (again, using the encoded values directly is no good because they'll be all over the place). Unique indices are contiguous again.
	std::map<funnelEncoding, funnelID> funnelTranslation;
	//vector giving the funnel ID for each individual
	std::vector<funnelID> lineFunnelIDs;
	//vector giving the encoded value for each individual
	std::vector<funnelEncoding> lineFunnelEncodings;
	//vector giving the encoded value for each value in allFunnels
	std::vector<funnelEncoding> allFunnelEncodings;
	funnelsToUniqueValues(funnelTranslation, lineFunnelIDs, lineFunnelEncodings, allFunnelEncodings, lineFunnels, allFunnels, nFounders);
	int nDifferentFunnels = lineFunnelEncodings.size();
	
	unsigned int maxAlleles = recoded.maxAlleles;
	if(maxAlleles > 64)
	{
		throw std::runtime_error("Internal error - Cannot have more than 64 alleles per marker");
	}

	markerPatternsToUniqueValuesArgs markerPatternData;
	markerPatternData.nFounders = nFounders;
	markerPatternData.nMarkers = nMarkers;
	markerPatternData.recodedFounders = recodedFounders;
	markerPatternData.recodedHetData = recodedHetData;
	markerPatternsToUniqueValues(markerPatternData);

	//Intermediate results. These give the most likely paths from the start of the chromosome to a marker, assuming some value for the underlying founder at the marker
	Rcpp::IntegerMatrix intermediate(nFounders, maxChromosomeMarkers);
	int cumulativeMarkerCounter = 0;

	xMajorMatrix<expandedProbabilitiesType> intercrossingHaplotypeProbabilities(maxChromosomeMarkers-1, maxAIGenerations, maxSelfing - minSelfing+1);
	rowMajorMatrix<expandedProbabilitiesType> funnelHaplotypeProbabilities(maxChromosomeMarkers-1, maxSelfing - minSelfing + 1);

	//We'll do a dispath based on whether or not we have infinite generations of selfing. Which requires partial template specialization, which requires a struct/class
	viterbiAlgorithm<nFounders, infiniteSelfing> viterbi(intercrossingHaplotypeProbabilities, funnelHaplotypeProbabilities, maxChromosomeMarkers);
	viterbi.recodedHetData = recodedHetData;
	viterbi.recodedFounders = recodedFounders;
	viterbi.recodedFinals = recodedFinals;
	viterbi.lineFunnelIDs = &lineFunnelIDs;
	viterbi.lineFunnelEncodings = &lineFunnelEncodings;
	viterbi.intercrossingGenerations = &intercrossingGenerations;
	viterbi.selfingGenerations = &selfingGenerations;
	viterbi.results = results;

	//Now actually run the Viterbi algorithm. To cut down on memory usage we run a single chromosome at a time
	for(int chromosomeCounter = 0; chromosomeCounter < map.size(); chromosomeCounter++)
	{
		//Generate haplotype probability data. 
		Rcpp::NumericVector positions = Rcpp::as<Rcpp::NumericVector>(map(chromosomeCounter));
		Rcpp::NumericVector recombinationFractions = haldaneToRf(diff(positions));
		for(int markerCounter = 0; markerCounter < recombinationFractions.size(); markerCounter++)
		{
			int markerPattern1 = markerPatternData.markerPatternIDs[markerCounter];
			int markerPattern2 = markerPatternData.markerPatternIDs[markerCounter + 1];
			markerData& markerPattern1Data = markerPatternData.allMarkerPatterns[markerPattern1];
			markerData& markerPattern2Data = markerPatternData.allMarkerPatterns[markerPattern2];
			double recombination = recombinationFractions(markerCounter);
			for(int selfingGenerationCounter = minSelfing; selfingGenerationCounter <= maxSelfing; selfingGenerationCounter++)
			{
				expandedGenotypeProbabilities<nFounders, infiniteSelfing>::noIntercross(funnelHaplotypeProbabilities(markerCounter, selfingGenerationCounter - minSelfing), recombination, selfingGenerationCounter, allFunnelEncodings.size());
			}
			for(int selfingGenerationCounter = minSelfing; selfingGenerationCounter <= maxSelfing; selfingGenerationCounter++)
			{
				for(int intercrossingGenerations = std::max(1, minAIGenerations); intercrossingGenerations <= maxAIGenerations; intercrossingGenerations++)
				{
					expandedGenotypeProbabilities<nFounders, infiniteSelfing>::withIntercross(intercrossingHaplotypeProbabilities(markerCounter, intercrossingGenerations - minAIGenerations, selfingGenerationCounter - minSelfing), intercrossingGenerations, recombination, selfingGenerationCounter, allFunnelEncodings.size());
				}
			}

		}
		//dispatch based on whether we have infinite generations of selfing or not. 
		viterbi.apply(cumulativeMarkerCounter, cumulativeMarkerCounter+positions.size());
		cumulativeMarkerCounter += positions.size();
	}
}
template<int nFounders> void imputedFoundersInternal1(Rcpp::IntegerMatrix founders, Rcpp::IntegerMatrix finals, Rcpp::S4 pedigree, Rcpp::List hetData, Rcpp::List map, Rcpp::IntegerMatrix results, bool infiniteSelfing)
{
	if(infiniteSelfing)
	{
		imputedFoundersInternal2<nFounders, true>(founders, finals, pedigree, hetData, map, results);
	}
	else
	{
		imputedFoundersInternal2<nFounders, false>(founders, finals, pedigree, hetData, map, results);
	}
}
SEXP imputeFounders(SEXP geneticData_sexp, SEXP map_sexp)
{
BEGIN_RCPP
	Rcpp::S4 geneticData;
	try
	{
		geneticData = Rcpp::as<Rcpp::S4>(geneticData_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData must be an S4 object of class geneticData");
	}

	Rcpp::IntegerMatrix founders;
	try
	{
		founders = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("founders"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@founders must be an integer matrix");
	}

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@finals must be an integer matrix");
	}

	Rcpp::S4 pedigree;
	try
	{
		pedigree = Rcpp::as<Rcpp::S4>(geneticData.slot("pedigree"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@pedigree must be an S4 object");
	}

	std::string pedigreeSelfingSlot;
	try
	{
		pedigreeSelfingSlot = Rcpp::as<std::string>(pedigree.slot("selfing"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@pedigree@selfing must be a string");
	}
	bool infiniteSelfing;
	if(pedigreeSelfingSlot == "infinite")
	{
		infiniteSelfing = true;
	}
	else if(pedigreeSelfingSlot == "auto")
	{
		infiniteSelfing = false;
	}
	else
	{
		throw std::runtime_error("Input geneticData@pedigree@selfing must be \"infinite\" or \"auto\"");
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

	Rcpp::List map;
	try
	{
		map = Rcpp::as<Rcpp::List>(map_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input map must be a list");
	}

	std::vector<std::string> foundersMarkers = Rcpp::as<std::vector<std::string> >(Rcpp::colnames(founders));
	std::vector<std::string> finalsMarkers = Rcpp::as<std::vector<std::string> >(Rcpp::colnames(finals));

	std::vector<std::string> mapMarkers;
	mapMarkers.reserve(foundersMarkers.size());
	int maxChromosomeMarkers = 0;
	for(int i = 0; i < map.size(); i++)
	{
		Rcpp::NumericVector chromosome;
		try
		{
			chromosome = Rcpp::as<Rcpp::NumericVector>(map(i));
		}
		catch(...)
		{
			throw std::runtime_error("Input map must be a list of numeric vectors");
		}
		Rcpp::CharacterVector chromosomeMarkers = chromosome.names();
		mapMarkers.insert(mapMarkers.end(), chromosomeMarkers.begin(), chromosomeMarkers.end());
		maxChromosomeMarkers = std::max(maxChromosomeMarkers, (int)chromosomeMarkers.size());
	}
	if(mapMarkers.size() != foundersMarkers.size() || !std::equal(mapMarkers.begin(), mapMarkers.end(), foundersMarkers.begin()))
	{
		throw std::runtime_error("Map was inconsistent with the markers in the geneticData object");
	}
	if(mapMarkers.size() != finalsMarkers.size() || !std::equal(mapMarkers.begin(), mapMarkers.end(), finalsMarkers.begin()))
	{
		throw std::runtime_error("Map was inconsistent with the markers in the geneticData object");
	}

	Rcpp::Function nFoundersFunc("nFounders");
	int nFounders = Rcpp::as<int>(nFoundersFunc(geneticData));
	int nFinals = finals.nrow();
	Rcpp::IntegerMatrix results(nFinals, mapMarkers.size());
	if(nFounders == 2)
	{
		imputedFoundersInternal1<2>(founders, finals, pedigree, hetData, map, results, infiniteSelfing);
	}
	else if(nFounders == 4)
	{
		imputedFoundersInternal1<4>(founders, finals, pedigree, hetData, map, results, infiniteSelfing);
	}
	else if(nFounders == 8)
	{
		imputedFoundersInternal1<8>(founders, finals, pedigree, hetData, map, results, infiniteSelfing);
	}
	else if(nFounders == 16)
	{
		imputedFoundersInternal1<16>(founders, finals, pedigree, hetData, map, results, infiniteSelfing);
	}
	else
	{
		throw std::runtime_error("Number of founders must be 2, 4, 8 or 16");
	}
	return results;
END_RCPP
}

