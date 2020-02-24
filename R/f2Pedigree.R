#' @title Generate an F2 pedigree which starts from inbred founders
#'
#' @description
#' Generate an F2 pedigree which starts from inbred founders
#' 
#' @param populationSize The size of the generated population. 
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
#' @examples 
#' pedigree <- f2Pedigree(1000)
#' #This pedigree is automatically marked as involving finite generations of selfing. 
#' selfing(pedigree)
f2Pedigree <- function(populationSize)
{
	pedigree <- twoParentPedigree(initialPopulationSize = 1, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = populationSize)
	pedigree@selfing <- "finite"
	return(pedigree)
}
