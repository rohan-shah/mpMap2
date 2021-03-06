% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roxygen.R
\docType{data}
\name{simulatedFourParentData}
\alias{simulatedFourParentData}
\alias{simulatedFourParentMap}
\title{Simulated data from a four-parent population.}
\description{
Simulated data from a four-parent population. Used in the examples given in the documentation.
}
\examples{
set.seed(1)
#This data was generated by the following script
pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000,
     selfingGenerations = 6, intercrossingGenerations = 0)
#Assume infinite generations of selfing in subsequent analysis
selfing(pedigree) <- "infinite"
#Generate random map
simulatedFourParentMap <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, 
	include.x = FALSE)
#Simulate data
simulatedFourParentData <- simulateMPCross(map = simulatedFourParentMap, pedigree = pedigree, 
	mapFunction = haldane, seed = 1L)
}
