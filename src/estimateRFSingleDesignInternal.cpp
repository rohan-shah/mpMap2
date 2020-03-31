#include "recodeFoundersFinalsHets.h"
#include "intercrossingAndSelfingGenerations.h"
#include "alleleDataErrors.h"
#include "estimateRFCheckFunnels.h"
#include "recodeHetsAsNA.h"
#include "warnings.h"
#include "matrixChunks.h"
#include "estimateRFSingleDesignInternal.h"
#include "constructLookupTable.h"
#include "intercrossingHaplotypeToMarker.h"
#include "funnelHaplotypeToMarker.h"
#include "probabilities.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#ifdef _OPENMP
#include "mpMap2_openmp.h"
#include <omp.h>
#endif
#include "matrix.h"
bool toInternalArgs(estimateRFSingleDesignArgs&& args, estimateRFSingleDesignInternalArgs& internal_args, std::string& error)
{
	error = "";
	std::stringstream ss;
	//work out the number of intercrossing generations
	int nFounders = args.founders.nrow(), nFinals = args.finals.nrow(), nMarkers = args.finals.ncol();
	std::vector<int> intercrossingGenerations, selfingGenerations;
	getIntercrossingAndSelfingGenerations(args.pedigree, args.finals, nFounders, intercrossingGenerations, selfingGenerations);
	bool hasAIC = nFinals > 0 && *std::max_element(intercrossingGenerations.begin(), intercrossingGenerations.end()) > 0;

	/*Check that all the observed marker values are potentially valid (ignoring pedigree). That is, is every observed value for the finals consistent with something in the hetData object?*/
	Rcpp::List codingErrors = listCodingErrors(args.founders, args.finals, args.hetData);
	std::vector<std::string> warnings, errors;
	codingErrorsToStrings(codingErrors, errors, args.finals, Rcpp::as<Rcpp::List>(args.hetData), 6);
	for(std::size_t errorIndex = 0; errorIndex < errors.size() && errorIndex < 6; errorIndex++)
	{
		ss << errors[errorIndex] << std::endl;
	}
	if(errors.size() > 0)
	{
		error = ss.str();
		return false;
	}

	//Check that everything has proper funnels - For the case of the lines which are just selfing, we just check that one funnel. For AIC lines, we check the funnels of all the parent lines
	//allFunnels stores a vector of all the funnels involved in any way. lineFunnels stores a value per line, specifying the funnel per line, if the line is not an intercrossing line. If it is a dummy value is inserted. 
	std::vector<funnelType> allFunnels, lineFunnels;
	{
		estimateRFCheckFunnels(args.finals, args.founders,  Rcpp::as<Rcpp::List>(args.hetData), args.pedigree, intercrossingGenerations, warnings, errors, allFunnels, lineFunnels);
		for(std::size_t errorIndex = 0; errorIndex < errors.size() && errorIndex < 6; errorIndex++)
		{
			ss << errors[errorIndex] << std::endl;;
		}
		if(errors.size() > 0)
		{
			error = ss.str();
			return false;
		}
		for(std::size_t warningIndex = 0; warningIndex < warnings.size() && warningIndex < 6; warningIndex++)
		{
			Rcpp::Rcout << warnings[warningIndex] << std::endl;
		}
		if(warnings.size() > 6)
		{
			Rcpp::Rcout << "Supressing further funnel warnings" << std::endl;
		}
	}
	//re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles
	//We do this to make it easier to identify markers with identical segregation patterns. recodedFounders = column major matrix
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	Rcpp::List recodedHetData(nMarkers);
	recodedHetData.attr("names") = args.hetData.attr("names");
	recodedFinals.attr("dimnames") = args.finals.attr("dimnames");
	recodedFounders.attr("dimnames") = args.founders.attr("dimnames");
	
	recodeDataStruct recoded;
	recoded.recodedFounders = recodedFounders;
	recoded.recodedFinals = recodedFinals;
	recoded.founders = args.founders;
	recoded.finals = args.finals;
	recoded.hetData = args.hetData;
	recoded.recodedHetData = recodedHetData;
	try
	{
		recodeFoundersFinalsHets(recoded);
	}
	catch(std::invalid_argument& argument)
	{
		throw std::runtime_error("Invalid input, please run validObject on the input mpcross object for more information");
	}

	//We need to assign a unique ID to each marker pattern - Where by pattern we mean the combination of hetData and founder alleles. Note that this is possible because we just recoded everything to a standardised format.
	//Marker IDs are guaranteed to be contiguous numbers starting from 0 - So the set of all valid [0, markerPatterns.size()]. 
	//Note that markerPatternID is defined in unitTypes.h. It's just an integer (and automatically convertible to an integer), but it's represented by a unique type - This stops us from confusing it with an ordinary integer.
	markerPatternsToUniqueValuesArgs markerPatternConversionArgs;
	markerPatternConversionArgs.nFounders = nFounders;
	markerPatternConversionArgs.nMarkers = nMarkers;
	markerPatternConversionArgs.recodedFounders = recodedFounders;
	markerPatternConversionArgs.recodedHetData = recodedHetData;
	markerPatternsToUniqueValues(markerPatternConversionArgs);

	unsigned int maxAlleles = 0;
	for(std::vector<int>::const_iterator i = internal_args.markerRows->begin(); i != internal_args.markerRows->end(); i++) 
	{
		markerPatternID currentMarkerID = markerPatternConversionArgs.markerPatternIDs[*i];
		maxAlleles = std::max(maxAlleles, (unsigned int)markerPatternConversionArgs.allMarkerPatterns[currentMarkerID].nObservedValues);
	}
	for(std::vector<int>::const_iterator i = internal_args.markerColumns->begin(); i != internal_args.markerColumns->end(); i++)
	{
		markerPatternID currentMarkerID = markerPatternConversionArgs.markerPatternIDs[*i];
		maxAlleles = std::max(maxAlleles, (unsigned int)markerPatternConversionArgs.allMarkerPatterns[currentMarkerID].nObservedValues);
	}
	if(maxAlleles > 64)
	{
		error = "To limit compilation time, cannot have more than 64 alleles per marker. Contact the author if you need this limit relaxed. ";
		return false;
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
	
	//In the case of infinite selfing, we've still allowed hets up to this point. But we want to ignore any potential hets in the final analysis. 
	bool infiniteSelfing = Rcpp::as<std::string>(args.pedigree.slot("selfing")) == "infinite";
	if(infiniteSelfing)
	{
		bool foundHets, foundHetEncodings;
		replaceHetsWithNA(recodedFounders, recodedFinals, recodedHetData, foundHets, foundHetEncodings);
		Rcpp::Function warning("warning");
		if(foundHets)
		{
			try
			{
				warning(hetWarning);
			}
			catch(...)
			{}
		}
		std::fill(selfingGenerations.begin(), selfingGenerations.end(), 0);
	}
	
	internal_args.finals = recodedFinals;
	internal_args.founders = recodedFounders;
	internal_args.pedigree = args.pedigree;
	internal_args.intercrossingGenerations.swap(intercrossingGenerations);
	internal_args.selfingGenerations.swap(selfingGenerations);
	internal_args.lineWeights.swap(args.lineWeights);
	internal_args.markerPatternData.swap(markerPatternConversionArgs);
	internal_args.hasAI = hasAIC;
	internal_args.lineFunnelIDs.swap(lineFunnelIDs);
	internal_args.lineFunnelEncodings.swap(lineFunnelEncodings);
	internal_args.allFunnelEncodings.swap(allFunnelEncodings);

	for(std::vector<int>::const_iterator i = internal_args.markerRows->begin(); i != internal_args.markerRows->end(); i++) internal_args.rowPatterns.push_back(internal_args.markerPatternData.markerPatternIDs[*i]);
	for(std::vector<int>::const_iterator i = internal_args.markerColumns->begin(); i != internal_args.markerColumns->end(); i++) internal_args.columnPatterns.push_back(internal_args.markerPatternData.markerPatternIDs[*i]);

	std::sort(internal_args.rowPatterns.begin(), internal_args.rowPatterns.end());
	std::sort(internal_args.columnPatterns.begin(), internal_args.columnPatterns.end());

	internal_args.rowPatterns.erase(std::unique(internal_args.rowPatterns.begin(), internal_args.rowPatterns.end()), internal_args.rowPatterns.end());
	internal_args.columnPatterns.erase(std::unique(internal_args.columnPatterns.begin(), internal_args.columnPatterns.end()), internal_args.columnPatterns.end());

	for(std::vector<int>::iterator rowPatternIterator = internal_args.rowPatterns.begin(); rowPatternIterator != internal_args.rowPatterns.end(); rowPatternIterator++)
	{
		for(std::vector<int>::iterator columnPatternIterator = internal_args.columnPatterns.begin(); columnPatternIterator != internal_args.columnPatterns.end(); columnPatternIterator++)
		{
			int i = *rowPatternIterator;
			int j = *columnPatternIterator;
			if(i > j) std::swap(i, j);
			{
				internal_args.uniquePatternPairs.push_back(std::make_pair(i, j));
			}
		}
	}
	std::sort(internal_args.uniquePatternPairs.begin(), internal_args.uniquePatternPairs.end());
	internal_args.uniquePatternPairs.erase(std::unique(internal_args.uniquePatternPairs.begin(), internal_args.uniquePatternPairs.end()), internal_args.uniquePatternPairs.end());
	return true;
}
template<int nFounders, bool infiniteSelfing> bool estimateRFSingleDesignInternal2(estimateRFSingleDesignInternalArgs& args)
{
	integerMatrix finals = args.finals;
	int nFinals = finals.nRow;

	int maxAIGenerations = *std::max_element(args.intercrossingGenerations.begin(), args.intercrossingGenerations.end());
	int minAIGenerations = getMinAIGenerations(&args.intercrossingGenerations);

	int maxSelfing = *std::max_element(args.selfingGenerations.begin(), args.selfingGenerations.end());
	int minSelfing = *std::min_element(args.selfingGenerations.begin(), args.selfingGenerations.end());

	const int maxAlleles = 64;
	
	std::size_t nRecombLevels = args.recombinationFractions.size(), nDifferentFunnels = args.lineFunnelEncodings.size();
	std::vector<double>::const_iterator halfIterator = std::find(args.recombinationFractions.begin(), args.recombinationFractions.end(), 0.5);
	int halfIndex = (int)std::distance(args.recombinationFractions.begin(), halfIterator);

	//In order to determine if a marker combination is informative, we use a much finer numerical grid.
	const int nFinerPoints = N_FINER_POINTS;
	std::vector<double> finerRecombLevels(nFinerPoints);
	for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
	{
		finerRecombLevels[recombCounter] = 0.5 * ((double)recombCounter) / ((double)nFinerPoints - 1.0);
	}

	typedef std::array<double, compressedProbabilities<nFounders, infiniteSelfing>::nDifferentProbs> compressedProbabilitiesType;
	rowMajorMatrix<compressedProbabilitiesType> finerFunnelHaplotypeProbabilities(nFinerPoints, maxSelfing-minSelfing+1);
	rowMajorMatrix<compressedProbabilitiesType> funnelHaplotypeProbabilities(nRecombLevels, maxSelfing - minSelfing + 1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int recombCounter = 0; recombCounter < (int)nRecombLevels; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(funnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), args.recombinationFractions[recombCounter], selfingGenerations, nDifferentFunnels);
		}
		for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
		{
			genotypeProbabilitiesNoIntercross<nFounders, infiniteSelfing>(finerFunnelHaplotypeProbabilities(recombCounter, selfingGenerations - minSelfing), finerRecombLevels[recombCounter], selfingGenerations, nDifferentFunnels);
		}
	}

	//Similarly for the intercrossing generation haplotype probabilities
	xMajorMatrix<compressedProbabilitiesType> intercrossingHaplotypeProbabilities(nRecombLevels, maxAIGenerations - minAIGenerations + 1, maxSelfing - minSelfing+1);
	xMajorMatrix<compressedProbabilitiesType> finerIntercrossingHaplotypeProbabilities(nFinerPoints, maxAIGenerations - minAIGenerations + 1, maxSelfing - minSelfing+1);
	for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
	{
		for(int aiCounter = minAIGenerations; aiCounter <= maxAIGenerations; aiCounter++)
		{
			for(int recombCounter = 0; recombCounter < (int)nRecombLevels; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(intercrossingHaplotypeProbabilities(recombCounter, aiCounter-minAIGenerations, selfingGenerations-minSelfing), aiCounter, args.recombinationFractions[recombCounter], selfingGenerations, nDifferentFunnels);
			}
			for(int recombCounter = 0; recombCounter < nFinerPoints; recombCounter++)
			{
				genotypeProbabilitiesWithIntercross<nFounders, infiniteSelfing>(finerIntercrossingHaplotypeProbabilities(recombCounter, aiCounter - minAIGenerations, selfingGenerations - minSelfing), aiCounter, finerRecombLevels[recombCounter], selfingGenerations, nDifferentFunnels);
			}
		}
	}
	const std::vector<markerPatternID>& markerPatternIDs = args.markerPatternData.markerPatternIDs;

	unsigned long long done = 0;
	const R_xlen_t product1 = maxAlleles*(maxSelfing-minSelfing + 1) *(nDifferentFunnels + maxAIGenerations - minAIGenerations+1);
	const R_xlen_t product2 = (maxSelfing - minSelfing + 1) *(nDifferentFunnels + maxAIGenerations - minAIGenerations + 1);
	const R_xlen_t product3 = nDifferentFunnels + maxAIGenerations - minAIGenerations + 1;
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		std::vector<array2<maxAlleles> > markerProbabilities(nFinerPoints);
		std::vector<int> table(maxAlleles*product1);
		std::vector<double> results(nRecombLevels);
		singleMarkerPairData<maxAlleles> thisMarkerPairData(nRecombLevels, nDifferentFunnels, std::max(maxAIGenerations - minAIGenerations+1, 0), maxSelfing - minSelfing + 1);
#ifdef _OPENMP
		#pragma omp for
#endif
		for(std::vector<std::pair<int, int> >::iterator patternPairIterator = args.uniquePatternPairs.begin(); patternPairIterator < args.uniquePatternPairs.end(); patternPairIterator++)
		{
			int markerPatternID1 = patternPairIterator->first;
			int markerPatternID2 = patternPairIterator->second;
				
			markerData& firstMarkerPatternData = args.markerPatternData.allMarkerPatterns[markerPatternID1];
			markerData& secondMarkerPatternData = args.markerPatternData.allMarkerPatterns[markerPatternID2];

			for(int selfingCounter = minSelfing; selfingCounter <= maxSelfing; selfingCounter++)
			{
				//Compute marker probabilities for a finer grid. If me seem to see a repeated probability model (numerically, up to a tolerance), then in that particular situtation this pair of markers is no good
				for(int funnelCounter = 0; funnelCounter < (int)nDifferentFunnels; funnelCounter++)
				{
					funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<false>(finerFunnelHaplotypeProbabilities, &(markerProbabilities[0]), args.lineFunnelEncodings[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing);
					thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues, finerRecombLevels);
				}
				for(int intercrossingGeneration = minAIGenerations; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
				{
					intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<false>(finerIntercrossingHaplotypeProbabilities, &(markerProbabilities[0]), intercrossingGeneration - minAIGenerations, firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing, args.allFunnelEncodings[0]);
					thisMarkerPairData.allowableAI(intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing) = isValid<maxAlleles>(markerProbabilities, nFinerPoints, firstMarkerPatternData.nObservedValues, secondMarkerPatternData.nObservedValues, finerRecombLevels);
				}
				//The next two loops relate to the input recombination fractions
				for(int intercrossingGeneration = minAIGenerations; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
				{
					array2<maxAlleles>* markerProbabilitiesThisIntercrossing = &(thisMarkerPairData.perAIGenerationData(0, intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing));
					if(thisMarkerPairData.allowableAI(intercrossingGeneration-minAIGenerations, selfingCounter - minSelfing))
					{
						intercrossingHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<true>(intercrossingHaplotypeProbabilities, markerProbabilitiesThisIntercrossing, intercrossingGeneration-minAIGenerations, firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing, args.allFunnelEncodings[0]);
					}
				}
				for(int funnelCounter = 0; funnelCounter < (int)nDifferentFunnels; funnelCounter++)
				{
					array2<maxAlleles>* markerProbabilitiesThisFunnel = &(thisMarkerPairData.perFunnelData(0, funnelCounter, selfingCounter - minSelfing));
					memset(markerProbabilitiesThisFunnel->values, 0, sizeof(array2<maxAlleles>::values));
					if(thisMarkerPairData.allowableFunnel(funnelCounter, selfingCounter - minSelfing))
					{
						funnelHaplotypeToMarker<nFounders, maxAlleles, infiniteSelfing>::template convert<true>(funnelHaplotypeProbabilities, markerProbabilitiesThisFunnel, args.lineFunnelEncodings[funnelCounter], firstMarkerPatternData, secondMarkerPatternData, selfingCounter - minSelfing);
					}
				}
			}
			std::function<bool(int)> predicate = [&markerPatternIDs,markerPatternID1,&markerPatternID2](int row)
			{
				return markerPatternIDs[row] == markerPatternID1 || markerPatternIDs[row] == markerPatternID2;
			};
			std::function<bool(int, int)> jointPredicate = [&markerPatternIDs,markerPatternID1,&markerPatternID2](int row, int column)
			{
				return (markerPatternIDs[row] == markerPatternID1 && markerPatternIDs[column] == markerPatternID2) || (markerPatternIDs[column] == markerPatternID1 && markerPatternIDs[row] == markerPatternID2);
			};
			triangularIteratorPredicates currentPosition(*args.markerRows, *args.markerColumns, predicate, jointPredicate);
			unsigned long long doneThisThread = 0;
			while(!currentPosition.isDone())
			{
				std::pair<int, int> markerIndices = currentPosition.get();
				int markerCounterRow = markerIndices.first, markerCounterColumn = markerIndices.second;
				
				int markerPatternIDRow = args.markerPatternData.markerPatternIDs[markerCounterRow];
				int markerPatternIDColumn = args.markerPatternData.markerPatternIDs[markerCounterColumn];
				bool swap = markerPatternIDRow > markerPatternIDColumn;
				if(swap) std::swap(markerPatternIDRow, markerPatternIDColumn);
				
#ifndef NDEBUG
				if(markerPatternID1 != markerPatternIDRow || markerPatternID2 != markerPatternIDColumn)
				{
					throw std::runtime_error("Internal error");
				}
#endif
				std::fill(table.begin(), table.end(), 0);
				for(int finalCounter = 0; finalCounter < (int)nFinals; finalCounter++)
				{
					int marker1Value = finals(finalCounter, markerCounterRow);
					int marker2Value = finals(finalCounter, markerCounterColumn);
					//If necessary swap the data
					if(swap) std::swap(marker1Value, marker2Value);
					if(marker1Value != NA_INTEGER && marker2Value != NA_INTEGER)
					{
						int intercrossingGenerations = args.intercrossingGenerations[finalCounter];
						int selfingGenerations = args.selfingGenerations[finalCounter];
						if(intercrossingGenerations == 0)
						{
							funnelID currentLineFunnelID = args.lineFunnelIDs[finalCounter];
							table[marker1Value*product1 + marker2Value*product2 + (selfingGenerations - minSelfing)*product3 + currentLineFunnelID]++;
						}
						else
						{
							table[marker1Value*product1 + marker2Value*product2 + (selfingGenerations - minSelfing)*product3 + nDifferentFunnels + intercrossingGenerations - minAIGenerations]++;
						}
					}
				}
				for(int recombCounter = 0; recombCounter < (int)nRecombLevels; recombCounter++)
				{
					double contribution = 0;
					for(int selfingGenerations = minSelfing; selfingGenerations <= maxSelfing; selfingGenerations++)
					{
						for(int marker1Value = 0; marker1Value < firstMarkerPatternData.nObservedValues; marker1Value++)
						{
							for(int marker2Value = 0; marker2Value < secondMarkerPatternData.nObservedValues; marker2Value++)
							{
								for(int intercrossingGenerations = minAIGenerations; intercrossingGenerations <= maxAIGenerations; intercrossingGenerations++)
								{
									int count = table[marker1Value*product1 + marker2Value * product2 + (selfingGenerations - minSelfing)*product3 + nDifferentFunnels + intercrossingGenerations - minAIGenerations];
									//This continue statement is important, because otherwise we may end up with 0 * -Inf, which results in -Inf
									if(count == 0) continue;
									bool allowable = thisMarkerPairData.allowableAI(intercrossingGenerations-minAIGenerations, selfingGenerations - minSelfing);
									if(allowable)
									{
										array2<maxAlleles>& perMarkerGenotypeValues = thisMarkerPairData.perAIGenerationData(recombCounter, intercrossingGenerations-minAIGenerations, selfingGenerations - minSelfing);
										contribution += count * perMarkerGenotypeValues.values[marker1Value][marker2Value];
									}
								}
								for(int funnelID = 0; funnelID < (int)nDifferentFunnels; funnelID++)
								{
									int count = table[marker1Value*product1 + marker2Value * product2 + (selfingGenerations - minSelfing)*product3 + funnelID];
									//This continue statement is important, because otherwise we may end up with 0 * -Inf, which results in -Inf
									if(count == 0) continue;
									bool allowable = thisMarkerPairData.allowableFunnel(funnelID, selfingGenerations - minSelfing);
									if(allowable)
									{
										array2<maxAlleles>& perMarkerGenotypeValues = thisMarkerPairData.perFunnelData(recombCounter, funnelID, selfingGenerations - minSelfing);
										contribution += count * perMarkerGenotypeValues.values[marker1Value][marker2Value];
									}
								}
							}
						}
					}
					if(contribution != contribution || contribution == -std::numeric_limits<double>::infinity()) results[recombCounter] = -std::numeric_limits<double>::infinity();
					else results[recombCounter] = contribution;
				}
				std::vector<double>::iterator maxIterator = std::max_element(results.begin(), results.end());
				double max = *maxIterator;
				double min = *std::min_element(results.begin(), results.end());
				double currentLod;
				int currentTheta = 0;
				if(max == 0 && min == 0)
				{
					max = currentLod = std::numeric_limits<double>::quiet_NaN();
					currentTheta = 0xff;
				}
				else
				{
					currentTheta = (int)std::distance(results.begin(), maxIterator);
					currentLod = max - results[halfIndex];
				}
				unsigned long long destinationCounter = currentPosition.getFlatIndex();
				args.theta[destinationCounter] = currentTheta;
				if(args.lod) args.lod[destinationCounter] = currentLod;
				if(args.lkhd) args.lkhd[destinationCounter] = max;
				doneThisThread++;
				currentPosition.next();
				destinationCounter++;
			}
#ifdef _OPENMP
			#pragma omp critical
#endif
			{
				done += doneThisThread;
#ifdef _OPENMP
				if(omp_get_thread_num() == 0)
#endif
				{
					args.updateProgress(done);
				}
			}
		}
	}
	return true;
}
template<bool infiniteSelfing> bool estimateRFSingleDesignInternal1(estimateRFSingleDesignInternalArgs& args)
{
	int nFounders = args.founders.nrow();
	if(nFounders == 2)
	{
		return estimateRFSingleDesignInternal2<2, infiniteSelfing>(args);
	}
	else if(nFounders == 4)
	{
		return estimateRFSingleDesignInternal2<4, infiniteSelfing>(args);
	}
	else if(nFounders == 8)
	{
		return estimateRFSingleDesignInternal2<8, infiniteSelfing>(args);
	}
	else if(nFounders == 16)
	{
		return estimateRFSingleDesignInternal2<16, infiniteSelfing>(args);
	}
	else
	{
		Rprintf("Number of founders must be 2, 4, 8 or 16\n");
		return false;
	}
	return true;
}
bool estimateRFSingleDesignInternal(estimateRFSingleDesignInternalArgs& args)
{
	bool infiniteSelfing = Rcpp::as<std::string>(args.pedigree.slot("selfing")) == "infinite";
	if(infiniteSelfing)
	{
		std::fill(args.selfingGenerations.begin(), args.selfingGenerations.end(), 0);
		return estimateRFSingleDesignInternal1<true>(args);
	}
	else return estimateRFSingleDesignInternal1<false>(args);
	
}
