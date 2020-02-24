setClass("removeHets", contains="NULL")
#' @export
#' @title Remove heterozygotes
#' @description Remove all heterozygotes from dataset
#' @details This function can be used to remove all heterozygotes from an \code{mpcross} object. Information about how pairs of different marker alleles are encoded as genotypes is discarded, and all observations of heterozygotes will be marked as \code{NA}. Any information calculated based on the genetic data (imputed IBD genotypes, IBD probabilities) will be discarded.
#' @return An object of internal class \code{removeHets}, which can be combined with an object of class \code{mpcross} using the addition operator. 
#' @examples
#' pedigree <- eightParentPedigreeImproperFunnels(initialPopulationSize = 10,
#'      selfingGenerations = 1, nSeeds = 1)
#' #Generate map
#' map <- qtl::sim.map()
#' #Simulate data
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
#' finals(cross)[1:5, 1:5]
#' hetData(cross)[[1]]
#' cross <- cross + removeHets()
#' finals(cross)[1:5, 1:5]
#' hetData(cross)[[1]]
removeHets <- function()
{
	return(new("removeHets"))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "removeHets"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Removing hets will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	for(i in 1:length(e1@geneticData))
	{
		newResults <- .Call("removeHets", e1@geneticData[[i]]@founders, e1@geneticData[[i]]@finals, e1@geneticData[[i]]@hetData, PACKAGE="mpMap2")
		e1@geneticData[[i]]@finals <- newResults$finals
		names(newResults$hetData) <- names(e1@geneticData[[i]]@hetData)
		e1@geneticData[[i]]@hetData <- new("hetData", newResults$hetData)
	}
	return(e1)
})
