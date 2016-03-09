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
		resultsMatrix <- .Call("imputeFounders", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, hetrozygoteMissingProb, PACKAGE="mpMap2")
		dimnames(resultsMatrix) <- dimnames(mpcrossMapped@geneticData[[i]]@finals)
		mpcrossMapped@geneticData[[i]]@imputed <- resultsMatrix
	}
	return(mpcrossMapped)
}
