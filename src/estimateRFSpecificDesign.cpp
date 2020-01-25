#include "estimateRFSpecificDesign.h"
#include <math.h>
#include "intercrossingAndSelfingGenerations.h"
#include "estimateRFCheckFunnels.h"
#include <map>
#include <set>
#include "recodeFoundersFinalsHets.h"
#include <memory>
#include "constructLookupTable.h"
#include "probabilities2.h"
#include "probabilities4.h"
#include "probabilities8.h"
#include "probabilities16.h"
#include "alleleDataErrors.h"
#include "recodeHetsAsNA.h"
#include "estimateRF.h"
#include "matrixChunks.h"
#ifdef USE_OPENMP
#include "mpMap2_openmp.h"
#include <omp.h>
#endif
#include "warnings.h"
#include "getMinAIGenerations.h"
template<int nFounders, int maxAlleles, bool infiniteSelfing> bool estimateRFSpecificDesign(rfhaps_internal_args& args, unsigned long long& progressCounter)
{
	std::size_t nFinals = args.finals.nrow(), nRecombLevels = args.recombinationFractions.size();
	if(nFinals == 0) return true;

	std::vector<double>& lineWeights = args.lineWeights;
	Rcpp::List finalDimNames = args.finals.attr("dimnames");
	Rcpp::CharacterVector finalNames = finalDimNames[0];

	int nMarkerPatternIDs = (int)args.markerPatternData.allMarkerPatterns.size();
	int minSelfing = *std::min_element(args.selfingGenerations.begin(), args.selfingGenerations.end());
	int minAIGenerations = getMinAIGenerations(&args.intercrossingGenerations);

	//This is basically just a huge lookup table
	allMarkerPairData<maxAlleles> computedContributions(nMarkerPatternIDs);
	
	constructLookupTableArgs<maxAlleles> lookupArgs(computedContributions, args.markerPatternData, args.rowPatterns, args.columnPatterns);
	lookupArgs.recombinationFractions = &args.recombinationFractions;
	lookupArgs.lineFunnelEncodings = &args.lineFunnelEncodings;
	lookupArgs.intercrossingGenerations = &args.intercrossingGenerations;
	lookupArgs.selfingGenerations = &args.selfingGenerations;
	lookupArgs.allFunnelEncodings = &args.allFunnelEncodings;
	constructLookupTable<nFounders, maxAlleles, infiniteSelfing>(lookupArgs);

	//We parallelise this array, even though it's over an iterator not an integer. So we use an integer and use that to work out how many steps forwards we need to move the iterator. We assume that the values are strictly increasing, otherwise this will never work. 
	//Use this to only call setTxtProgressBar every 10 calls to updateProgress. Probably no point in updating status more frequently than that.
	unsigned long long updateProgressCounter = 0;
#ifdef USE_OPENMP
	#pragma omp parallel 
#endif
	{
		triangularIterator indexIterator = args.startPosition;
		unsigned long long previousCounter = 0;
#ifdef USE_OPENMP
		#pragma omp for schedule(dynamic)
#endif
		for(unsigned long long counter = 0; counter < args.valuesToEstimateInChunk; counter++)
		{
			signed long long difference = counter - previousCounter;
			if(difference < 0LL) throw std::runtime_error("Internal error");
			while(difference > 0LL) 
			{
				indexIterator.next();
				difference--;
			}
			previousCounter = counter;

			std::pair<int, int> markerIndices = indexIterator.get();
			int markerCounterRow = markerIndices.first, markerCounterColumn = markerIndices.second;

			int markerPatternID1 = args.markerPatternData.markerPatternIDs[markerCounterRow];
			int markerPatternID2 = args.markerPatternData.markerPatternIDs[markerCounterColumn];

			singleMarkerPairData<maxAlleles>& markerPairData = *(computedContributions(markerPatternID1, markerPatternID2));
			//We only calculated tabels for markerPattern1 <= markerPattern2. So if we want things the other way around we have to swap the data for markers 1 and 2 later on. 
			bool swap = markerPatternID1 > markerPatternID2;
			for(int recombCounter = 0; recombCounter < (int)nRecombLevels; recombCounter++)
			{
				for(int finalCounter = 0; finalCounter < (int)nFinals; finalCounter++)
				{
					int marker1Value = args.finals(finalCounter, markerCounterRow);
					int marker2Value = args.finals(finalCounter, markerCounterColumn);
					//If necessary swap the data
					if(swap) std::swap(marker1Value, marker2Value);
					if(marker1Value != NA_INTEGER && marker2Value != NA_INTEGER)
					{
						double contribution = 0;
						bool allowable = false;
						int intercrossingGenerations = args.intercrossingGenerations[finalCounter];
						int selfingGenerations = args.selfingGenerations[finalCounter];
						if(intercrossingGenerations == 0)
						{
							funnelID currentLineFunnelID = args.lineFunnelIDs[finalCounter];
							allowable = markerPairData.allowableFunnel(currentLineFunnelID, selfingGenerations - minSelfing);
							if(allowable)
							{
								array2<maxAlleles>& perMarkerGenotypeValues = markerPairData.perFunnelData(recombCounter, currentLineFunnelID, selfingGenerations - minSelfing);
								contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
							}
						}
						else if(intercrossingGenerations > 0)
						{
							allowable = markerPairData.allowableAI(intercrossingGenerations-minAIGenerations, selfingGenerations - minSelfing);
							if(allowable)
							{
								array2<maxAlleles>& perMarkerGenotypeValues = markerPairData.perAIGenerationData(recombCounter, intercrossingGenerations-minAIGenerations, selfingGenerations-minSelfing);
								contribution = perMarkerGenotypeValues.values[marker1Value][marker2Value];
							}
						}
						//We get an NA from trying to take the logarithm of zero - That is, this parameter is completely impossible for the given data, so put in -Inf
						if(contribution != contribution || contribution == -std::numeric_limits<double>::infinity()) args.result[(long)counter *(long)nRecombLevels + (long)recombCounter] = -std::numeric_limits<double>::infinity();
						else if(contribution != 0 && allowable) args.result[(long)counter * (long)nRecombLevels + (long)recombCounter] += lineWeights[finalCounter] * contribution;
					}
				}
			}
#ifdef USE_OPENMP
			#pragma omp critical
#endif
			{
				progressCounter++;
			}
#ifdef USE_OPENMP
			if(omp_get_thread_num() == 0)
#endif
			{
				updateProgressCounter++;
				if(updateProgressCounter % 100 == 0) args.updateProgress(progressCounter);
			}
		}
	}
	return true;
}
template<int nFounders, int maxAlleles, bool infiniteSelfing> bool estimateRFSpecificDesignNoLineWeights(rfhaps_internal_args& args, unsigned long long& progressCounter)
{
	std::size_t nFinals = args.finals.nrow(), nRecombLevels = args.recombinationFractions.size();
	if(nFinals == 0) return true;

	std::size_t nDifferentFunnels = args.lineFunnelEncodings.size();
	Rcpp::List finalDimNames = args.finals.attr("dimnames");
	Rcpp::CharacterVector finalNames = finalDimNames[0];

	int nMarkerPatternIDs = (int)args.markerPatternData.allMarkerPatterns.size();
	int maxAIGenerations = *std::max_element(args.intercrossingGenerations.begin(), args.intercrossingGenerations.end());
	int minAIGenerations = getMinAIGenerations(&args.intercrossingGenerations);
	int minSelfing = *std::min_element(args.selfingGenerations.begin(), args.selfingGenerations.end());
	int maxSelfing = *std::max_element(args.selfingGenerations.begin(), args.selfingGenerations.end());

	//This is basically just a huge lookup table
	allMarkerPairData<maxAlleles> computedContributions(nMarkerPatternIDs);
	
	constructLookupTableArgs<maxAlleles> lookupArgs(computedContributions, args.markerPatternData, args.rowPatterns, args.columnPatterns);
	lookupArgs.recombinationFractions = &args.recombinationFractions;
	lookupArgs.lineFunnelEncodings = &args.lineFunnelEncodings;
	lookupArgs.intercrossingGenerations = &args.intercrossingGenerations;
	lookupArgs.selfingGenerations = &args.selfingGenerations;
	lookupArgs.allFunnelEncodings = &args.allFunnelEncodings;
	constructLookupTable<nFounders, maxAlleles, infiniteSelfing>(lookupArgs);

	const R_xlen_t product1 = maxAlleles*(maxSelfing-minSelfing + 1) *(nDifferentFunnels + maxAIGenerations - minAIGenerations+1);
	const R_xlen_t product2 = (maxSelfing - minSelfing + 1) *(nDifferentFunnels + maxAIGenerations - minAIGenerations + 1);
	const R_xlen_t product3 = nDifferentFunnels + maxAIGenerations - minAIGenerations + 1;

	//We parallelise this array, even though it's over an iterator not an integer. So we use an integer and use that to work out how many steps forwards we need to move the iterator. We assume that the values are strictly increasing, otherwise this will never work.
	//Use this to only call setTxtProgressBar every 10 calls to updateProgress. Probably no point in updating status more frequently than that.
	unsigned long long updateProgressCounter = 0;
#ifdef USE_OPENMP
	#pragma omp parallel 
#endif
	{
		triangularIterator indexIterator = args.startPosition;
		//Indexing is of the form table[allele1 * product1 + allele2*product2 + selfingGenerations * product3 + (ai OR funnel)]. Funnels come first. 
		std::vector<int> table(maxAlleles*product1);

		int previousCounter = 0;
#ifdef USE_OPENMP
		#pragma omp for schedule(dynamic)
#endif
		for(int counter = 0; counter < (int)args.valuesToEstimateInChunk; counter++)
		{
			std::fill(table.begin(), table.end(), 0);
			int difference = counter - previousCounter;
			if(difference < 0) throw std::runtime_error("Internal error");
			while(difference > 0) 
			{
				indexIterator.next();
				difference--;
			}
			previousCounter = counter;

			std::pair<int, int> markerIndices = indexIterator.get();
			int markerCounterRow = markerIndices.first, markerCounterColumn = markerIndices.second;

			int markerPatternID1 = args.markerPatternData.markerPatternIDs[markerCounterRow];
			int markerPatternID2 = args.markerPatternData.markerPatternIDs[markerCounterColumn];

			bool swap = markerPatternID1 > markerPatternID2;

			int markerPatternID1Copied = markerPatternID1, markerPatternID2Copied = markerPatternID2;
			if(swap) std::swap(markerPatternID1Copied, markerPatternID2Copied);

			singleMarkerPairData<maxAlleles>& markerPairData = *(computedContributions(markerPatternID1Copied, markerPatternID2Copied));
			//We only calculated tabels for markerPattern1 <= markerPattern2. So if we want things the other way around we have to swap the data for markers 1 and 2 later on. 
			for(int finalCounter = 0; finalCounter < (int)nFinals; finalCounter++)
			{
				int marker1Value = args.finals(finalCounter, markerCounterRow);
				int marker2Value = args.finals(finalCounter, markerCounterColumn);
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
					else if(intercrossingGenerations > 0)
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
					for(int marker1Value = 0; marker1Value < maxAlleles; marker1Value++)
					{
						for(int marker2Value = 0; marker2Value < maxAlleles; marker2Value++)
						{
							for(int intercrossingGenerations = minAIGenerations; intercrossingGenerations <= maxAIGenerations; intercrossingGenerations++)
							{
								int count = table[marker1Value*product1 + marker2Value * product2 + (selfingGenerations - minSelfing)*product3 + nDifferentFunnels + intercrossingGenerations - minAIGenerations];
								//This continue statement is important, because otherwise we may end up with 0 * -Inf, which results in -Inf
								if(count == 0) continue;
								bool allowable = markerPairData.allowableAI(intercrossingGenerations-minAIGenerations, selfingGenerations - minSelfing);
								if(allowable)
								{
									array2<maxAlleles>& perMarkerGenotypeValues = markerPairData.perAIGenerationData(recombCounter, intercrossingGenerations-minAIGenerations, selfingGenerations - minSelfing);
									contribution += count * perMarkerGenotypeValues.values[marker1Value][marker2Value];
								}
							}
							for(int funnelID = 0; funnelID < (int)nDifferentFunnels; funnelID++)
							{
								int count = table[marker1Value*product1 + marker2Value * product2 + (selfingGenerations - minSelfing)*product3 + funnelID];
								//This continue statement is important, because otherwise we may end up with 0 * -Inf, which results in -Inf
								if(count == 0) continue;
								bool allowable = markerPairData.allowableFunnel(funnelID, selfingGenerations - minSelfing);
								if(allowable)
								{
									array2<maxAlleles>& perMarkerGenotypeValues = markerPairData.perFunnelData(recombCounter, funnelID, selfingGenerations - minSelfing);
									contribution += count * perMarkerGenotypeValues.values[marker1Value][marker2Value];
								}
							}

						}
					}
				}
				//We get an NA from trying to take the logarithm of zero - That is, this parameter is completely impossible for the given data, so put in -Inf
				if(contribution != contribution || contribution == -std::numeric_limits<double>::infinity()) args.result[(long)counter *(long)nRecombLevels + (long)recombCounter] = -std::numeric_limits<double>::infinity();
				else args.result[(long)counter * (long)nRecombLevels + (long)recombCounter] += contribution;
			}
#ifdef USE_OPENMP
			#pragma omp critical
#endif
			{
				progressCounter++;
			}
#ifdef USE_OPENMP
			if(omp_get_thread_num() == 0)
#endif
			{
				updateProgressCounter++;
				if(updateProgressCounter % 100 == 0) args.updateProgress(progressCounter);
			}
		}
	}
	return true;
}
template<int nFounders, int maxAlleles, bool infiniteSelfing> bool estimateRFSpecificDesign3(rfhaps_internal_args& args, unsigned long long& counter)
{
	for(std::vector<double>::iterator i = args.lineWeights.begin(); i != args.lineWeights.end(); i++)
	{
		if(*i != 1) return estimateRFSpecificDesign<nFounders, maxAlleles, infiniteSelfing>(args, counter);
	}
	return estimateRFSpecificDesignNoLineWeights<nFounders, maxAlleles, infiniteSelfing>(args, counter);
}
template<int nFounders, int maxAlleles> bool estimateRFSpecificDesignInternal2(rfhaps_internal_args& args, unsigned long long& counter)
{
	bool infiniteSelfing = Rcpp::as<std::string>(args.pedigree.slot("selfing")) == "infinite";
	if(infiniteSelfing)
	{
		std::fill(args.selfingGenerations.begin(), args.selfingGenerations.end(), 0);
		return estimateRFSpecificDesign3<nFounders, maxAlleles, true>(args, counter);
	}
	else return estimateRFSpecificDesign3<nFounders, maxAlleles, false>(args, counter);
}
//here we transfer maxAlleles over to the templated parameter section - This can make a BIG difference to memory usage if this is smaller, and it's going into a type so it has to be templated.
template<int nFounders> bool estimateRFSpecificDesignInternal1(rfhaps_internal_args& args, unsigned long long& counter)
{
	//for i in `seq 1 64`; do echo -e "case $i:\n\t\treturn estimateRFSpecificDesignInternal2<nFounders, $i>(args, counter);"; done
	switch(args.maxAlleles)
	{
		case 1:
		case 2:
			return estimateRFSpecificDesignInternal2<nFounders, 2>(args, counter);
		case 3:
			return estimateRFSpecificDesignInternal2<nFounders, 3>(args, counter);
		case 4:
			return estimateRFSpecificDesignInternal2<nFounders, 4>(args, counter);
		case 5:
		case 6:
			return estimateRFSpecificDesignInternal2<nFounders, 6>(args, counter);
		case 7:
		case 8:
			return estimateRFSpecificDesignInternal2<nFounders, 8>(args, counter);
		case 9:
		case 10:
			return estimateRFSpecificDesignInternal2<nFounders, 10>(args, counter);
		case 11:
		case 12:
			return estimateRFSpecificDesignInternal2<nFounders, 12>(args, counter);
		case 13:
		case 14:
			return estimateRFSpecificDesignInternal2<nFounders, 14>(args, counter);
		case 15:
		case 16:
			return estimateRFSpecificDesignInternal2<nFounders, 16>(args, counter);
		case 17:
		case 18:
			return estimateRFSpecificDesignInternal2<nFounders, 18>(args, counter);
		case 19:
		case 20:
			return estimateRFSpecificDesignInternal2<nFounders, 20>(args, counter);
		case 21:
		case 22:
			return estimateRFSpecificDesignInternal2<nFounders, 22>(args, counter);
		case 23:
		case 24:
			return estimateRFSpecificDesignInternal2<nFounders, 24>(args, counter);
		case 25:
		case 26:
			return estimateRFSpecificDesignInternal2<nFounders, 26>(args, counter);
		case 27:
		case 28:
			return estimateRFSpecificDesignInternal2<nFounders, 28>(args, counter);
		case 29:
		case 30:
			return estimateRFSpecificDesignInternal2<nFounders, 30>(args, counter);
		case 31:
		case 32:
			return estimateRFSpecificDesignInternal2<nFounders, 32>(args, counter);
		case 33:
		case 34:
			return estimateRFSpecificDesignInternal2<nFounders, 34>(args, counter);
		case 35:
		case 36:
			return estimateRFSpecificDesignInternal2<nFounders, 36>(args, counter);
		case 37:
		case 38:
			return estimateRFSpecificDesignInternal2<nFounders, 38>(args, counter);
		case 39:
		case 40:
			return estimateRFSpecificDesignInternal2<nFounders, 40>(args, counter);
		case 41:
		case 42:
			return estimateRFSpecificDesignInternal2<nFounders, 42>(args, counter);
		case 43:
		case 44:
			return estimateRFSpecificDesignInternal2<nFounders, 44>(args, counter);
		case 45:
		case 46:
			return estimateRFSpecificDesignInternal2<nFounders, 46>(args, counter);
		case 47:
		case 48:
			return estimateRFSpecificDesignInternal2<nFounders, 48>(args, counter);
		case 49:
		case 50:
			return estimateRFSpecificDesignInternal2<nFounders, 50>(args, counter);
		case 51:
		case 52:
			return estimateRFSpecificDesignInternal2<nFounders, 52>(args, counter);
		case 53:
		case 54:
			return estimateRFSpecificDesignInternal2<nFounders, 54>(args, counter);
		case 55:
		case 56:
			return estimateRFSpecificDesignInternal2<nFounders, 56>(args, counter);
		case 57:
		case 58:
			return estimateRFSpecificDesignInternal2<nFounders, 58>(args, counter);
		case 59:
		case 60:
			return estimateRFSpecificDesignInternal2<nFounders, 60>(args, counter);
		case 61:
		case 62:
			return estimateRFSpecificDesignInternal2<nFounders, 62>(args, counter);
		case 63:
		case 64:
			return estimateRFSpecificDesignInternal2<nFounders, 64>(args, counter);
		default:
			throw std::runtime_error("Internal error");
	}
}
unsigned long long estimateLookup(rfhaps_internal_args& internal_args)
{
	if(internal_args.finals.nrow() == 0)
	{
		return 0ULL;
	}
	int maxAIGenerations = *std::max_element(internal_args.intercrossingGenerations.begin(), internal_args.intercrossingGenerations.end());
	int minAIGenerations = getMinAIGenerations(&internal_args.intercrossingGenerations);
	int minSelfing = *std::min_element(internal_args.selfingGenerations.begin(), internal_args.selfingGenerations.end());
	int maxSelfing = *std::max_element(internal_args.selfingGenerations.begin(), internal_args.selfingGenerations.end());
	std::size_t nDifferentFunnels = internal_args.lineFunnelEncodings.size();
	std::size_t nRecombLevels = internal_args.recombinationFractions.size();
	
	std::size_t arraySize;
	//for i in `seq 1 64`; do echo -e "\t\tcase $i:\n\t\t\tarraySize = sizeof(array2<$i>);\n\t\t\tbreak;"; done > tmp
	switch(internal_args.maxAlleles)
	{
		case 1:
		case 2:
			arraySize = sizeof(array2<2>);
			break;
		case 3:
			arraySize = sizeof(array2<3>);
			break;
		case 4:
			arraySize = sizeof(array2<4>);
			break;
		case 5:
		case 6:
			arraySize = sizeof(array2<6>);
			break;
		case 7:
		case 8:
			arraySize = sizeof(array2<8>);
			break;
		case 9:
		case 10:
			arraySize = sizeof(array2<10>);
			break;
		case 11:
		case 12:
			arraySize = sizeof(array2<12>);
			break;
		case 13:
		case 14:
			arraySize = sizeof(array2<14>);
			break;
		case 15:
		case 16:
			arraySize = sizeof(array2<16>);
			break;
		case 17:
		case 18:
			arraySize = sizeof(array2<18>);
			break;
		case 19:
		case 20:
			arraySize = sizeof(array2<20>);
			break;
		case 21:
		case 22:
			arraySize = sizeof(array2<22>);
			break;
		case 23:
		case 24:
			arraySize = sizeof(array2<24>);
			break;
		case 25:
		case 26:
			arraySize = sizeof(array2<26>);
			break;
		case 27:
		case 28:
			arraySize = sizeof(array2<28>);
			break;
		case 29:
		case 30:
			arraySize = sizeof(array2<30>);
			break;
		case 31:
		case 32:
			arraySize = sizeof(array2<32>);
			break;
		case 33:
		case 34:
			arraySize = sizeof(array2<34>);
			break;
		case 35:
		case 36:
			arraySize = sizeof(array2<36>);
			break;
		case 37:
		case 38:
			arraySize = sizeof(array2<38>);
			break;
		case 39:
		case 40:
			arraySize = sizeof(array2<40>);
			break;
		case 41:
		case 42:
			arraySize = sizeof(array2<42>);
			break;
		case 43:
		case 44:
			arraySize = sizeof(array2<44>);
			break;
		case 45:
		case 46:
			arraySize = sizeof(array2<46>);
			break;
		case 47:
		case 48:
			arraySize = sizeof(array2<48>);
			break;
		case 49:
		case 50:
			arraySize = sizeof(array2<50>);
			break;
		case 51:
		case 52:
			arraySize = sizeof(array2<52>);
			break;
		case 53:
		case 54:
			arraySize = sizeof(array2<54>);
			break;
		case 55:
		case 56:
			arraySize = sizeof(array2<56>);
			break;
		case 57:
		case 58:
			arraySize = sizeof(array2<58>);
			break;
		case 59:
		case 60:
			arraySize = sizeof(array2<60>);
			break;
		case 61:
		case 62:
			arraySize = sizeof(array2<62>);
			break;
		case 63:
		case 64:
			arraySize = sizeof(array2<64>);
			break;
		default:
			throw std::runtime_error("Internal error");
	}
	int nMarkerPairs = 0;
	for(std::vector<int>::iterator markerRowPattern = internal_args.rowPatterns.begin(); markerRowPattern != internal_args.rowPatterns.end(); markerRowPattern++)
	{
		std::vector<int>::iterator firstBiggerOrEqual = std::lower_bound(internal_args.columnPatterns.begin(), internal_args.columnPatterns.end(), *markerRowPattern);
		nMarkerPairs += (int)std::distance(firstBiggerOrEqual, internal_args.columnPatterns.end());
	}
	return nRecombLevels * (maxSelfing - minSelfing + 1) * (nDifferentFunnels + maxAIGenerations - minAIGenerations + 1) * arraySize * nMarkerPairs;
}
bool toInternalArgs(estimateRFSpecificDesignArgs&& args, rfhaps_internal_args& internal_args, std::string& error)
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
	
	//We don't use recoded.maxAlleles, because we want to restrict to the markers in markerRows and markerColumns
	unsigned int maxAlleles = 0; //= recoded.maxAlleles;
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
	//vector giving the encoded value for each funnelID in lineFunnelIDs
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
	internal_args.maxAlleles = maxAlleles;
	internal_args.lineFunnelIDs.swap(lineFunnelIDs);
	internal_args.lineFunnelEncodings.swap(lineFunnelEncodings);
	internal_args.allFunnelEncodings.swap(allFunnelEncodings);

	for(std::vector<int>::const_iterator i = internal_args.markerRows->begin(); i != internal_args.markerRows->end(); i++) internal_args.rowPatterns.push_back(internal_args.markerPatternData.markerPatternIDs[*i]);
	for(std::vector<int>::const_iterator i = internal_args.markerColumns->begin(); i != internal_args.markerColumns->end(); i++) internal_args.columnPatterns.push_back(internal_args.markerPatternData.markerPatternIDs[*i]);

	std::sort(internal_args.rowPatterns.begin(), internal_args.rowPatterns.end());
	std::sort(internal_args.columnPatterns.begin(), internal_args.columnPatterns.end());

	internal_args.rowPatterns.erase(std::unique(internal_args.rowPatterns.begin(), internal_args.rowPatterns.end()), internal_args.rowPatterns.end());
	internal_args.columnPatterns.erase(std::unique(internal_args.columnPatterns.begin(), internal_args.columnPatterns.end()), internal_args.columnPatterns.end());
	return true;
}
bool estimateRFSpecificDesign(rfhaps_internal_args& internal_args, unsigned long long& counter)
{
	int nFounders = internal_args.founders.nrow();
	if(nFounders == 2)
	{
		return estimateRFSpecificDesignInternal1<2>(internal_args, counter);
	}
	else if(nFounders == 4)
	{
		return estimateRFSpecificDesignInternal1<4>(internal_args, counter);
	}
	else if(nFounders == 8)
	{
		return estimateRFSpecificDesignInternal1<8>(internal_args, counter);
	}
	else if(nFounders == 16)
	{
		return estimateRFSpecificDesignInternal1<16>(internal_args, counter);
	}
	else
	{
		Rprintf("Number of founders must be 2, 4, 8 or 16\n");
		return false;
	}
	return true;
}
