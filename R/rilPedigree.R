#' @title Generate a two-parent RIL pedigree which starts from inbred founders
#' 
#' @description 
#' Generate a two-parent RIL pedigree which starts from inbred founders
#'
#' @param populationSize The size of the generated population
#' @export
rilPedigree <- function(populationSize, selfingGenerations)
{
	return(twoParentPedigree(initialPopulationSize = populationSize, selfingGenerations = selfingGenerations, intercrossingGenerations = 0, nSeeds = 1))
}
