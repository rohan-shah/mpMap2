#' @title Get positions of genetic markers
#' @description Get positions of genetic markers, on their respective chromosomes
#' @details Get positions of genetic markers in cM, on their respective chromosomes
#' @examples
#' map <- qtl::sim.map()
#' pedigree <- f2Pedigree(1000)
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
#' mappedCross <- mpcrossMapped(cross = cross, map = map)
#' getPositions(mappedCross, c("D13M3", "DXM1", "DXM3"))
#' @param mpcrossMapped The \code{mpcross} object containing the map of interest
#' @param markers The markers for which to get the positions 
#' @return A named vector of numbers, with names corresponding to the selected genetic markers, and values corresponding to genetic positions. 
#' @export
getPositions <- function(mpcrossMapped, markers)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(missing(markers) || !is.character(markers) || any(is.na(markers)))
	{
		stop("Input markers must be a character vector, without missing values")
	}
	markerPositions <- unlist(sapply(markers, function(marker)
	{
		na.omit(unlist(lapply(mpcrossMapped@map, function(x) x[marker])))
	}))
	names(markerPositions) <- markers
	return(markerPositions)
}
