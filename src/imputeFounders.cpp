#include "imputeFounders.h"
#include "intercrossingAndSelfingGenerations.h"
#include "recodeFoundersFinalsHets.h"
#include "matrices.h"
#include "probabilities.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#include "funnelsToUniqueValues.h"
#include "estimateRFCheckFunnels.h"
#include "markerPatternsToUniqueValues.h"
#include "intercrossingHaplotypeToMarker.h"
#include "funnelHaplotypeToMarker.h"
#include "viterbi.h"
#include "recodeHetsAsNA.h"
#include "warnings.h"
#include "joinMapWithExtra.h"
#include "haldaneToRf.h"
#include "generateKeys.h"
#include "getMinAIGenerations.h"
template<int nFounders, bool infiniteSelfing> void imputedFoundersInternal2(Rcpp::IntegerMatrix founders, Rcpp::IntegerMatrix finals, Rcpp::S4 pedigree, Rcpp::List hetData, Rcpp::IntegerMatrix results, Rcpp::IntegerMatrix resultsErrors, double homozygoteMissingProb, double heterozygoteMissingProb, double errorProb, Rcpp::IntegerMatrix key, positionData& allPositions, bool showProgress)
{
	//Work out maximum number of positions per chromosome
	int maxChromosomePositions = 0;
	for(std::size_t i = 0; i < allPositions.chromosomes.size(); i++)
	{
		positionData::chromosomeDescriptor& currentChromosome = allPositions.chromosomes[i];
		maxChromosomePositions = std::max(currentChromosome.end - currentChromosome.start, maxChromosomePositions);
	}

	typedef typename expandedProbabilities<nFounders, infiniteSelfing>::type expandedProbabilitiesType;
	//expandedProbabilitiesType haplotypeProbabilities;

	//Get out generations of selfing and intercrossing
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(pedigree, finals, nFounders, intercrossingGenerations, selfingGenerations);

	int maxSelfing = *std::max_element(selfingGenerations.begin(), selfingGenerations.end());
	int minSelfing = *std::min_element(selfingGenerations.begin(), selfingGenerations.end());
	int maxAIGenerations = *std::max_element(intercrossingGenerations.begin(), intercrossingGenerations.end());
	int minAIGenerations = getMinAIGenerations(&intercrossingGenerations);
	int nFinals = finals.nrow();
	int nMarkers = founders.ncol();

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
	try
	{
		recodeFoundersFinalsHets(recoded);
	}
	catch(std::invalid_argument& argument)
	{
		throw std::runtime_error("Invalid input, please run validObject on the input mpcross object for more information");
	}

	if(infiniteSelfing)
	{
		bool foundHets, foundHetEncodings;
		replaceHetsWithNA(recodedFounders, recodedFinals, recodedHetData, foundHets, foundHetEncodings);
		Rcpp::Function warning("warning");
		if(foundHets)
		{
			//Technically a warning could lead to an error if options(warn=2). This would be bad because it would break out of our code. This solution generates a c++ exception in that case, which we can then ignore. 
			try
			{
				warning(hetWarning);
			}
			catch(...)
			{}
		}
		std::fill(selfingGenerations.begin(), selfingGenerations.end(), 0);
	}

	//Get out the number of unique funnels. This is only needed because in the case of one funnel we assume a single funnel design and in the case of multiple funnels we assume a random funnels design
	std::vector<std::string> errors, warnings;
	std::vector<funnelType> allFunnels, lineFunnels;
	{
		estimateRFCheckFunnels(recodedFinals, recodedFounders, recodedHetData, pedigree, intercrossingGenerations, warnings, errors, allFunnels, lineFunnels);
		if(errors.size() > 0)
		{
			std::stringstream ss;
			for(std::size_t i = 0; i < errors.size(); i++)
			{
				ss << errors[i] << std::endl;
			}
			throw std::runtime_error(ss.str().c_str());
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
	for(std::vector<funnelEncoding>::iterator i = allFunnelEncodings.begin(); i != allFunnelEncodings.end(); i++)
	{
		funnelEncoding enc = *i;
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		std::sort(funnel, funnel + nFounders);
		if(std::unique(funnel, funnel + nFounders) != funnel + nFounders)
		{
			throw std::runtime_error("Duplicate founders in a funnel. The code for this case is not yet written");
		}
	}
	
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
	Rcpp::IntegerMatrix intermediate(nFounders, maxChromosomePositions);

	xMajorMatrix<expandedProbabilitiesType> logIntercrossingHaplotypeProbabilities(maxChromosomePositions, maxAIGenerations - minAIGenerations + 1, maxSelfing - minSelfing+1);
	rowMajorMatrix<expandedProbabilitiesType> logFunnelHaplotypeProbabilities(maxChromosomePositions, maxSelfing - minSelfing + 1);

	//The single loci probabilities are different depending on whether there are zero or one generations of intercrossing. But once you have non-zero generations, it doesn't matter how many
	std::vector<array2<nFounders> > logIntercrossingSingleLociHaplotypeProbabilities(maxSelfing - minSelfing+1), intercrossingSingleLociHaplotypeProbabilities(maxSelfing - minSelfing+1);
	std::vector<array2<nFounders> > logFunnelSingleLociHaplotypeProbabilities(maxSelfing - minSelfing + 1), funnelSingleLociHaplotypeProbabilities(maxSelfing - minSelfing + 1);

	int nFunnels = (int)allFunnelEncodings.size();
	//Generate single loci genetic data
	for(int selfingGenerationCounter = minSelfing; selfingGenerationCounter <= maxSelfing; selfingGenerationCounter++)
	{
		array2<nFounders>& funnelArray = funnelSingleLociHaplotypeProbabilities[selfingGenerationCounter - minSelfing];
		array2<nFounders>& intercrossingArray = intercrossingSingleLociHaplotypeProbabilities[selfingGenerationCounter - minSelfing];
		singleLocusGenotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(funnelArray, selfingGenerationCounter, nFunnels);
		singleLocusGenotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(intercrossingArray, selfingGenerationCounter, nFunnels);

		array2<nFounders>& logFunnelArray = logFunnelSingleLociHaplotypeProbabilities[selfingGenerationCounter - minSelfing];
		array2<nFounders>& logIntercrossingArray = logIntercrossingSingleLociHaplotypeProbabilities[selfingGenerationCounter - minSelfing];
		//Take logarithms
		for(int i = 0; i < nFounders; i++)
		{
			for(int j = 0; j < nFounders; j++)
			{
				if(funnelArray.values[i][j] == 0) logFunnelArray.values[i][j] = -std::numeric_limits<double>::infinity();
				else logFunnelArray.values[i][j] = log(funnelArray.values[i][j]);
				if(intercrossingArray.values[i][j] == 0) logIntercrossingArray.values[i][j] = -std::numeric_limits<double>::infinity();
				else logIntercrossingArray.values[i][j] = log(intercrossingArray.values[i][j]);
			}
		}
	}

	//We'll do a dispath based on whether or not we have infinite generations of selfing. Which requires partial template specialization, which requires a struct/class
	viterbiAlgorithm<nFounders, infiniteSelfing> viterbi(markerPatternData, maxChromosomePositions, allPositions);
	viterbi.logIntercrossingHaplotypeProbabilities = &logIntercrossingHaplotypeProbabilities;
	viterbi.logFunnelHaplotypeProbabilities = &logFunnelHaplotypeProbabilities;
	viterbi.recodedHetData = recodedHetData;
	viterbi.recodedFounders = recodedFounders;
	viterbi.recodedFinals = recodedFinals;
	viterbi.lineFunnelIDs = &lineFunnelIDs;
	viterbi.lineFunnelEncodings = &lineFunnelEncodings;
	viterbi.intercrossingGenerations = &intercrossingGenerations;
	viterbi.selfingGenerations = &selfingGenerations;
	viterbi.results = results;
	viterbi.resultsErrors = resultsErrors;
	viterbi.key = key;
	viterbi.maxAlleles = maxAlleles; 
	viterbi.homozygoteMissingProb = homozygoteMissingProb;
	viterbi.heterozygoteMissingProb = heterozygoteMissingProb;
	viterbi.errorProb = errorProb;
	viterbi.logIntercrossingSingleLociHaplotypeProbabilities = &logIntercrossingSingleLociHaplotypeProbabilities;
	viterbi.logFunnelSingleLociHaplotypeProbabilities = &logFunnelSingleLociHaplotypeProbabilities;
	viterbi.intercrossingSingleLociHaplotypeProbabilities = &intercrossingSingleLociHaplotypeProbabilities;
	viterbi.funnelSingleLociHaplotypeProbabilities = &funnelSingleLociHaplotypeProbabilities;
	viterbi.minAIGenerations = minAIGenerations;
	viterbi.maxAIGenerations = maxAIGenerations;

	Rcpp::Function txtProgressBar("txtProgressBar");
	Rcpp::Function setTxtProgressBar("setTxtProgressBar");
	Rcpp::Function close("close");
	Rcpp::RObject barHandle;
	std::function<void(unsigned long long)> updateProgress = [](unsigned long long){};
	if(showProgress)
	{
		updateProgress = [&barHandle,nFinals,setTxtProgressBar](unsigned long long done)
		{
			try
			{
#ifdef CUSTOM_STATIC_RCPP
				setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)nFinals));
#else
				setTxtProgressBar(barHandle, (int)((double)(1000*done) / (double)nFinals));
#endif
			}
			catch(...)
			{
			}
		};
	}
	viterbi.updateProgress = updateProgress;

	//Now actually run the Viterbi algorithm. To cut down on memory usage we run a single chromosome at a time
	for(int chromosomeCounter = 0; chromosomeCounter < (int)allPositions.chromosomes.size(); chromosomeCounter++)
	{
		positionData::chromosomeDescriptor& currentChromosome = allPositions.chromosomes[chromosomeCounter];
		if(showProgress)
		{
			Rcpp::Rcout << "Begining chromosome " << currentChromosome.name << std::endl;
		}
		//Take differences
		std::vector<double> differences(currentChromosome.end - currentChromosome.start - 1);
		for(std::size_t i = 1; i < (std::size_t)(currentChromosome.end - currentChromosome.start); i++) differences[i-1] = allPositions.positions[i + currentChromosome.start] - allPositions.positions[i + currentChromosome.start -1];

		//Convert to recombination fractions
		std::vector<double> recombinationFractions;
		haldaneToRf(differences, recombinationFractions);
		//Generate haplotype probability data. 
		for(std::size_t markerCounter = 0; markerCounter < recombinationFractions.size(); markerCounter++)
		{
			double recombination = recombinationFractions[markerCounter];
			for(int selfingGenerationCounter = minSelfing; selfingGenerationCounter <= maxSelfing; selfingGenerationCounter++)
			{
				expandedGenotypeProbabilitiesUnphased<nFounders, infiniteSelfing, true>::noIntercross(logFunnelHaplotypeProbabilities(markerCounter, selfingGenerationCounter - minSelfing), recombination, selfingGenerationCounter, nFunnels);
			}
			for(int selfingGenerationCounter = minSelfing; selfingGenerationCounter <= maxSelfing; selfingGenerationCounter++)
			{
				for(int intercrossingGenerations =  minAIGenerations; intercrossingGenerations <= maxAIGenerations; intercrossingGenerations++)
				{
					expandedGenotypeProbabilitiesUnphased<nFounders, infiniteSelfing, true>::withIntercross(logIntercrossingHaplotypeProbabilities(markerCounter, intercrossingGenerations - minAIGenerations, selfingGenerationCounter - minSelfing), intercrossingGenerations, recombination, selfingGenerationCounter, nFunnels);
				}
			}

		}
		if(showProgress)
		{
			barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
		}
		//dispatch based on whether we have infinite generations of selfing or not. 
		viterbi.apply(currentChromosome.start, currentChromosome.end);
		if(showProgress)
		{
			close(barHandle);
		}
	}
	Rcpp::rownames(results) = Rcpp::rownames(finals);
	Rcpp::colnames(results) = Rcpp::wrap(allPositions.names);
}
template<int nFounders> void imputedFoundersInternal1(Rcpp::IntegerMatrix founders, Rcpp::IntegerMatrix finals, Rcpp::S4 pedigree, Rcpp::List hetData, Rcpp::IntegerMatrix results, Rcpp::IntegerMatrix resultsErrors, bool infiniteSelfing, double homozygoteMissingProb, double heterozygoteMissingProb, double errorProb, Rcpp::IntegerMatrix key, positionData& allPositions, bool showProgress)
{
	if(infiniteSelfing)
	{
		imputedFoundersInternal2<nFounders, true>(founders, finals, pedigree, hetData, results, resultsErrors, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
	}
	else
	{
		imputedFoundersInternal2<nFounders, false>(founders, finals, pedigree, hetData, results, resultsErrors, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
	}
}
SEXP imputeFounders(SEXP geneticData_sexp, SEXP map_sexp, SEXP homozygoteMissingProb_sexp, SEXP heterozygoteMissingProb_sexp, SEXP errorProb_sexp, SEXP extraPositions_sexp, SEXP showProgress_sexp)
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

	bool showProgress;
	try
	{
		showProgress = Rcpp::as<bool>(showProgress_sexp);
	}
	catch(...)
	{
		throw  std::runtime_error("Input showProgress must be a boolean");
	}
	//Check that row and column names are not null
	if(Rcpp::as<Rcpp::RObject>(Rcpp::rownames(founders)).isNULL())
	{
		throw std::runtime_error("Input founders must have row names");
	}
	if(Rcpp::as<Rcpp::RObject>(Rcpp::colnames(founders)).isNULL())
	{
		throw std::runtime_error("Input founders must have column names");
	}
	//Check that there are no values of NA
	{
		Rcpp::CharacterVector foundersMarkers = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(founders));
		for(std::size_t i = 0; i < (std::size_t)foundersMarkers.size(); i++)
		{
			if(Rcpp::CharacterVector::is_na(foundersMarkers[i])) throw std::runtime_error("Input founders cannot have a column named NA");
		}
		Rcpp::CharacterVector lineNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::rownames(founders));
		for(std::size_t i = 0; i < (std::size_t)lineNames.size(); i++)
		{
			if(Rcpp::CharacterVector::is_na(lineNames[i])) throw std::runtime_error("Input founders cannot have a row named NA");
		}
	}
	//Finally actually convert to a C++ type. 
	std::vector<std::string> foundersMarkers = Rcpp::as<std::vector<std::string> >(Rcpp::colnames(founders));

	Rcpp::IntegerMatrix finals;
	try
	{
		finals = Rcpp::as<Rcpp::IntegerMatrix>(geneticData.slot("finals"));
	}
	catch(...)
	{
		throw std::runtime_error("Input geneticData@finals must be an integer matrix");
	}
	//Check that row and column names are not null
	if(Rcpp::as<Rcpp::RObject>(Rcpp::rownames(finals)).isNULL())
	{
		throw std::runtime_error("Input finals must have row names");
	}
	if(Rcpp::as<Rcpp::RObject>(Rcpp::colnames(finals)).isNULL())
	{
		throw std::runtime_error("Input finals must have column names");
	}
	//Check that there are no values of NA
	{
		Rcpp::CharacterVector finalsMarkers = Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(finals));
		for(std::size_t i = 0; i < foundersMarkers.size(); i++)
		{
			if(Rcpp::CharacterVector::is_na(finalsMarkers[i])) throw std::runtime_error("Input finals cannot have a column named NA");
		}
		Rcpp::CharacterVector lineNames = Rcpp::as<Rcpp::CharacterVector>(Rcpp::rownames(finals));
		for(std::size_t i = 0; i < (std::size_t)lineNames.size(); i++)
		{
			if(Rcpp::CharacterVector::is_na(lineNames[i])) throw std::runtime_error("Input finals cannot have a row named NA");
		}
	}
	//Finally actually convert to a C++ type. 
	std::vector<std::string> finalsMarkers = Rcpp::as<std::vector<std::string> >(Rcpp::colnames(finals));
	std::vector<std::string> lineNames = Rcpp::as<std::vector<std::string> >(Rcpp::rownames(finals));

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
	else if(pedigreeSelfingSlot == "finite")
	{
		infiniteSelfing = false;
	}
	else
	{
		throw std::runtime_error("Input geneticData@pedigree@selfing must be \"infinite\" or \"finite\"");
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

	double homozygoteMissingProb;
	try
	{
		homozygoteMissingProb = Rcpp::as<double>(homozygoteMissingProb_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input homozygoteMissingProb must be a number between 0 and 1");
	}
	if(homozygoteMissingProb < 0 || homozygoteMissingProb > 1) throw std::runtime_error("Input homozygoteMissingProb must be a number between 0 and 1");

	double heterozygoteMissingProb;
	try
	{
		heterozygoteMissingProb = Rcpp::as<double>(heterozygoteMissingProb_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input heterozygoteMissingProb must be a number between 0 and 1");
	}
	if(heterozygoteMissingProb < 0 || heterozygoteMissingProb > 1) throw std::runtime_error("Input heterozygoteMissingProb must be a number between 0 and 1");

	double errorProb;
	try
	{
		errorProb = Rcpp::as<double>(errorProb_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input errorProb must be a number between 0 and 1");
	}
	if(errorProb < 0 || errorProb > 1) throw std::runtime_error("Input errorProb must be a number between 0 and 1");

	Rcpp::List extraPositions;
	try
	{
		extraPositions = Rcpp::as<Rcpp::List>(extraPositions_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input extraPositions must be a list");
	}

	//Unify the 

	Rcpp::Function nFoundersFunc("nFounders");
	int nFounders = Rcpp::as<int>(nFoundersFunc(geneticData));

	//Construct the key that takes pairs of founder values and turns them into encodings
	//We also want a version closer to the hetData format
	Rcpp::IntegerMatrix key, outputKey;
	generateKeys(key, outputKey, nFounders, infiniteSelfing);

	std::vector<std::string> mapMarkers;
	mapMarkers.reserve(foundersMarkers.size());
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
		std::vector<std::string> chromosomeMarkers = Rcpp::as<std::vector<std::string> >(chromosome.names());
		mapMarkers.insert(mapMarkers.end(), chromosomeMarkers.begin(), chromosomeMarkers.end());
	}
	if(mapMarkers.size() != foundersMarkers.size() || !std::equal(mapMarkers.begin(), mapMarkers.end(), foundersMarkers.begin()))
	{
		throw std::runtime_error("Map was inconsistent with the markers in the geneticData object");
	}
	if(mapMarkers.size() != finalsMarkers.size() || !std::equal(mapMarkers.begin(), mapMarkers.end(), finalsMarkers.begin()))
	{
		throw std::runtime_error("Map was inconsistent with the markers in the geneticData object");
	}

	positionData allPositions;
	joinMapWithExtra(map, extraPositions, allPositions);

	int nFinals = finals.nrow();
	Rcpp::IntegerMatrix results(nFinals, (int)allPositions.names.size());
	Rcpp::IntegerMatrix resultsErrors(nFinals, (int)allPositions.names.size());
	try
	{
		if(nFounders == 2)
		{
			imputedFoundersInternal1<2>(founders, finals, pedigree, hetData, results, resultsErrors, infiniteSelfing, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
		}
		else if(nFounders == 4)
		{
			imputedFoundersInternal1<4>(founders, finals, pedigree, hetData, results, resultsErrors, infiniteSelfing, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
		}
		else if(nFounders == 8)
		{
			imputedFoundersInternal1<8>(founders, finals, pedigree, hetData, results, resultsErrors, infiniteSelfing, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
		}
		else if(nFounders == 16)
		{
			imputedFoundersInternal1<16>(founders, finals, pedigree, hetData, results, resultsErrors, infiniteSelfing, homozygoteMissingProb, heterozygoteMissingProb, errorProb, key, allPositions, showProgress);
		}
		else
		{
			throw std::runtime_error("Number of founders must be 2, 4, 8 or 16");
		}
	}
	catch(impossibleDataException& err)
	{
		std::stringstream ss;
		ss << "Impossible data may have been detected for markers " << mapMarkers[err.marker] << " and " << mapMarkers[err.marker+1] << " for line " << lineNames[err.line] << ". Are these markers at the same location, and if so does this line have a recombination event between these markers?"; 
		throw std::runtime_error(ss.str().c_str());
	}
	return Rcpp::List::create(Rcpp::Named("data") = results, Rcpp::Named("key") = outputKey, Rcpp::Named("map") = allPositions.makeUnifiedMap(), Rcpp::Named("errors") = resultsErrors);
END_RCPP
}

