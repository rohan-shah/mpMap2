#' @export
computeGenotypeProbabilities <- function(mpcrossMapped, homozygoteMissingProb = 1, heterozygoteMissingProb = 1)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(homozygoteMissingProb < 0 || homozygoteMissingProb > 1)
	{
		stop("Input homozygoteMissingProb must be a value between 0 and 1")
	}
	if(heterozygoteMissingProb < 0 || heterozygoteMissingProb > 1)
	{
		stop("Input heterozygoteMissingProb must be a value between 0 and 1")
	}
	for(i in 1:length(mpcrossMapped@geneticData))
	{
		results <- .Call("computeGenotypeProbabilities", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, heterozygoteMissingProb, PACKAGE="mpMap2")
		resultsMatrix <- results$data
		nAlleles <- nrow(resultsMatrix) / nrow(mpcrossMapped@geneticData[[i]]@finals)
		colnames(resultsMatrix) <- colnames(mpcrossMapped@geneticData[[i]]@finals)
		rownames(resultsMatrix) <- unlist(lapply(rownames(mpcrossMapped@geneticData[[i]]@finals), function(lineName) paste0(lineName, " - ", 1:nAlleles)))
		mpcrossMapped@geneticData[[i]]@probabilities <- new("probabilities", data = resultsMatrix, key = results$key)
	}
	return(mpcrossMapped)
}
