#' @export
getIntercrossingAndSelfingGenerations <- function(cross)
{
	if(length(cross@geneticData) == 1)
	{
		result <- .Call("getIntercrossingAndSelfingGenerationsExport", cross@geneticData[[1]]@pedigree, cross@geneticData[[1]]@finals, PACKAGE="mpMap2")
		asMatrix <- do.call(cbind, result)
		rownames(asMatrix) <- rownames(finals(cross))
		return(asMatrix)
	}
	else
	{
		result <- lapply(cross@geneticData, 
			function(x)
			{
				asMatrix <- do.call(cbind, .Call("getIntercrossingAndSelfingGenerationsExport", x@pedigree, x@finals, PACKAGE="mpMap2"))
				rownames(asMatrix) <- rownames(x@finals)
				return(asMatrix)
			})
		return(result)
	}
}
