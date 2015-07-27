#' @title Special two-parent designs
#' @name specialTwoParent
#'
#' @description
#' Generate common two-parent designs
NULL
#' @param populationSize Final population size
#' @param selfingGenerations Number of generations of selfing to perform
#' @describeIn specialTwoParent Generate a pedigree for a RIL design
#' @export
rilPedigree <- function(populationSize, selfingGenerations)
{
	nonNegativeIntegerArgument(populationSize)
  	nonNegativeIntegerArgument(selfingGenerations)
	pedigree <- twoParentPedigree(initialPopulationSize = populationSize, selfingGenerations=selfingGenerations, nSeeds=1, intercrossingGenerations=0)
	return(pedigree)
}
#' @describeIn specialTwoParent Generate a pedigree for an F2 design
#' @param populationSize Final population size
#' @export
f2Pedigree <- function(populationSize)
{
	pedigree <- twoParentPedigree(initialPopulationSize = 1, selfingGenerations=1, nSeeds=populationSize, intercrossingGenerations=0)
	pedigree@selfing <- "auto"
	return(pedigree)	
}