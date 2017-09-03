#' @export
computeAllEpistaticChiSquared <- function(mpcrossMapped, showProgress = TRUE)
{
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Function computeAllEpistaticChiSquared can only be applied to a single design at a time")
	}
	if(is.null(mpcrossMapped@geneticData[[1]]@probabilities))
	{
		stop("Function computeAllEpistaticChiSquared requires probabilities to have been calculated")
	}
	epistatic <- .Call("computeAllEpistaticChiSquared", mpcrossMapped@geneticData[[1]]@probabilities, nFounders(mpcrossMapped), selfing(mpcrossMapped@geneticData[[1]]@pedigree) == "infinite", showProgress, PACKAGE="mpMap2")
	probabilitiesMap <- mpcrossMapped@geneticData[[1]]@probabilities@map
	rownames(epistatic) <- colnames(epistatic) <- unlist(lapply(probabilitiesMap, names))
	return(epistatic)
}
