#' @export
imputeFounders <- function(mpcrossMapped, homozygoteMissingProb = 1, hetrozygoteMissingProb = 1)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(homozygoteMissingProb < 0 || homozygoteMissingProb > 1)
	{
		stop("Input homozygoteMissingProb must be a value between 0 and 1")
	}
	if(hetrozygoteMissingProb < 0 || hetrozygoteMissingProb > 1)
	{
		stop("Input hetrozygoteMissingProb must be a value between 0 and 1")
	}
	for(i in 1:length(mpcrossMapped@geneticData))
	{
		results <- .Call("imputeFounders", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, hetrozygoteMissingProb, PACKAGE="mpMap2")
		resultsMatrix <- results$data
		dimnames(resultsMatrix) <- dimnames(mpcrossMapped@geneticData[[i]]@finals)
		mpcrossMapped@geneticData[[i]]@imputed <- new("imputed", data = resultsMatrix, key = results$key)
	}
	return(mpcrossMapped)
}
