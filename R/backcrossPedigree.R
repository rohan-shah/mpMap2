#' @title Generate a backcross pedigree which starts from inbred founders
#'
#' @description Generate a backcross pedigree which starts from inbred founders
#' @details Generate a backcross pedigree which starts from inbred founders
#' 
#' @param populationSize The size of the generated population. 
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @examples 
#' pedigree <- backcrossPedigree(1000)
#' #This pedigree is automatically marked as involving finite generations of selfing. 
#' selfing(pedigree)
#' @export
backcrossPedigree <- function(populationSize)
{
	entries <- 6 + populationSize
	father <- mother <- vector(mode = "integer", entries)
	observed <- vector(mode = "logical", entries)
	observed[7:(6+populationSize)] <- TRUE
	lineNames <- paste0("L", 1:entries)

	mother[1:4] <- father[1:4] <- 0L
	mother[5] <- father[5] <- 2L
	mother[6] <- 1L
	father[6] <- 2L
	mother[7:(6+populationSize)] <- 6L
	father[7:(6+populationSize)] <- 5L
	return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1:4L, observed = observed, selfing = "finite", warnImproperFunnels = FALSE))
}
