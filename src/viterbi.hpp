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
#include <limits>
template<int nFounders, bool infiniteSelfing> struct viterbiAlgorithm;
template<int nFounders> struct viterbiAlgorithm<nFounders, false>
{
	typedef typename expandedProbabilities<nFounders, false>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	rowMajorMatrix<double> intermediate1, intermediate2;
	Rcpp::IntegerMatrix results;
	std::vector<double> pathLengths1, pathLengths2;
	std::vector<double> working;
	xMajorMatrix<array2<16> >& intercrossingMarkerProbabilities;
	xMajorMatrix<array2<16> >& funnelMarkerProbabilities;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minSelfingGenerations;
	int maxSelfingGenerations;
	viterbiAlgorithm(xMajorMatrix<array2<16> >& intercrossingMarkerProbabilities, xMajorMatrix<array2<16> >& funnelMarkerProbabilities, int maxChromosomeSize)
		: intermediate1(nFounders, maxChromosomeSize), intermediate2(nFounders, maxChromosomeSize), pathLengths1(nFounders), pathLengths2(nFounders), working(nFounders), intercrossingMarkerProbabilities(intercrossingMarkerProbabilities), funnelMarkerProbabilities(funnelMarkerProbabilities)
	{}
	void apply(int start, int end)
	{
	}
};
template<int nFounders> struct viterbiAlgorithm<nFounders, true>
{
	typedef typename expandedProbabilities<nFounders, true>::type expandedProbabilitiesType;
	Rcpp::List recodedHetData;
	Rcpp::IntegerMatrix recodedFounders, recodedFinals;
	rowMajorMatrix<double> intermediate1, intermediate2;
	Rcpp::IntegerMatrix results;
	std::vector<double> pathLengths1, pathLengths2;
	std::vector<double> working;
	xMajorMatrix<array2<16> >& intercrossingMarkerProbabilities;
	xMajorMatrix<array2<16> >& funnelMarkerProbabilities;
	std::vector<funnelID>* lineFunnelIDs;
	std::vector<funnelEncoding>* lineFunnelEncodings;
	std::vector<int>* intercrossingGenerations;
	std::vector<int>* selfingGenerations;
	int minSelfingGenerations;
	int maxSelfingGenerations;
	viterbiAlgorithm(xMajorMatrix<array2<16> >& intercrossingMarkerProbabilities, xMajorMatrix<array2<16> >& funnelMarkerProbabilities, int maxChromosomeSize)
		: intermediate1(nFounders, maxChromosomeSize), intermediate2(nFounders, maxChromosomeSize), pathLengths1(nFounders), pathLengths2(nFounders), working(nFounders), intercrossingMarkerProbabilities(intercrossingMarkerProbabilities), funnelMarkerProbabilities(funnelMarkerProbabilities)
	{}
	void apply(int start, int end)
	{
		minSelfingGenerations = *std::min_element(selfingGenerations->begin(), selfingGenerations->end());
		maxSelfingGenerations = *std::max_element(selfingGenerations->begin(), selfingGenerations->end());
		int nFinals = recodedFinals.nrow(), nMarkers = recodedFinals.ncol();
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if((*intercrossingGenerations)[finalCounter] == 0)
			{
				applyFunnel(start, end, finalCounter, (*lineFunnelIDs)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			else
			{
				applyIntercrossing(start, end, finalCounter, (*intercrossingGenerations)[finalCounter], (*selfingGenerations)[finalCounter]);
			}
			std::vector<double>::iterator shortestPath = std::min_element(pathLengths1.begin(), pathLengths1.end());
			int shortestIndex = std::distance(pathLengths1.begin(), shortestPath);
			for(int i = 0; i < end - start; i++)
			{
				results(finalCounter, i) = intermediate1(shortestIndex, i);
			}
		}
	}
	void applyFunnel(int start, int end, int finalCounter, int funnelID, int selfingGenerations)
	{
		//Initialise the algorithm. For infinite generations of selfing, we don't need to bother with the hetData object, as there are no hets
		int markerValue = recodedFinals(finalCounter, start);
		funnelEncoding enc = (*lineFunnelEncodings)[finalCounter];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			intermediate1(founderCounter, 0) = founderCounter;
			pathLengths1[0] = 0;
			if(recodedFounders(founderCounter, start) != markerValue)
			{
				pathLengths1[nFounders] = -std::numeric_limits<double>::infinity();
			}
		}
		for(int markerCounter = start; markerCounter < end - 1; markerCounter++)
		{
			int previousMarkerValue = recodedFinals(finalCounter, markerCounter);
			markerValue = recodedFinals(finalCounter, markerCounter+1);
			//The founder at the next marker
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(funnel[founderCounter], markerCounter+1) == markerValue)
				{
					//Founder at the previous marker. 
					std::fill(working.begin(), working.end(), std::numeric_limits<double>::infinity());
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						if(recodedFounders(funnel[founderCounter2], markerCounter) == previousMarkerValue)
						{
							working[funnel[founderCounter2]] = pathLengths1[funnel[founderCounter2]] + funnelMarkerProbabilities(markerCounter, funnelID, selfingGenerations - minSelfingGenerations).values[founderCounter2][founderCounter];
						}
					}
					//Get the shortest one, and check that it's not negative infinity.
					std::vector<double>::iterator shortest = std::min_element(working.begin(), working.end());
					if(*shortest == std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
					int bestPrevious = (int)std::distance(working.begin(), shortest);
					
					memcpy(&(intermediate2(founderCounter, 0)), &(intermediate1(bestPrevious, 0)), sizeof(double)*(markerCounter - start + 1));
					intermediate2(founderCounter, markerCounter-start+1) = founderCounter;
					pathLengths2[founderCounter] = *shortest;
				}
				else
				{
					pathLengths2[founderCounter] = std::numeric_limits<double>::infinity();
				}
			}
			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
		}
	}
	void applyIntercrossing(int start, int end, int finalCounter, int intercrossingGenerations, int selfingGenerations)
	{
	}
};

