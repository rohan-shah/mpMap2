#' @title Generate an F2 pedigree which starts from inbred founders
#'
#' @description
#' Generate an F2 pedigree which starts from inbred founders
#' 
#' @param populationSize The size of the generated population
#' @export
f2Pedigree <- function(populationSize)
{
	return(twoParentPedigree(initialPopulationSize = populationSize, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1))
}
