#' @title Identify number of generations of intercrossing and selfing, per genetic line
#' @description Identify number of generations of intercrossing and selfing, per genetic line
#' @details Many structured populations consist of a number of generations of mixing, followed by a number of generations of intercrossing, followed by inbreeding. This function identifies the number of generations of selfing and intercrossing, for each genetic line, in the case of 4-way, 8-way or 16-way multi-parent design. 
#' @param cross The \code{mpcross} object containing the pedigree to be analysed. 
#' @return An integer matrix with two columns, giving the number of generations of selfing and intercrossing, for each genetic line. Or in the case of multiple experiments contained within a single object, a list of such matrices. 
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
