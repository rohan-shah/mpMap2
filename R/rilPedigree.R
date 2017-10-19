#' @title Generate a two-parent RIL pedigree which starts from inbred founders
#' 
#' @description 
#' Generate a two-parent RIL pedigree which starts from inbred founders
#'
#' @param populationSize The size of the generated population
#' @param selfingGenerations Number of generations of selfing. Specifying one generation leads to an F2 design. 
#' @export
rilPedigree <- function(populationSize, selfingGenerations)
{
	return(twoParentPedigree(initialPopulationSize = populationSize, selfingGenerations = selfingGenerations, intercrossingGenerations = 0, nSeeds = 1))
}
