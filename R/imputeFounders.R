#' @export
imputeFounders <- function(mpcrossMapped)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	for(i in 1:length(mpcrossMapped@geneticData))
	{
		resultsMatrix <- .Call("imputeFounders", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, PACKAGE="mpMap2")
		dimnames(resultsMatrix) <- dimnames(mpcrossMapped@geneticData[[i]]@finals)
		mpcrossMapped@geneticData[[i]]@imputed <- resultsMatrix
	}
	return(mpcrossMapped)
}
