#' @title Generate an F2 pedigree which starts from inbred founders
#'
#' @description
#' Generate an F2 pedigree which starts from inbred founders
#' 
#' @param populationSize The size of the generated population
#' @export
f2Pedigree <- function(populationSize)
{
	pedigree <- twoParentPedigree(initialPopulationSize = 1, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = populationSize)
	pedigree@selfing <- "auto"
	return(pedigree)
}
