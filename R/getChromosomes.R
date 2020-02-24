#' @title Get chromosome assignment per marker
#' @description Get chromosome assignment per marker from an \code{mpcross} object.
#' @details Extract a character vector, with names corresponding to markers, and values corresponding to the chromosome on which the named marker is located.
#' @param markers The markers for which we want the chromosomes
#' @param mpcrossMapped The object containing the map of interest
#' @return A character vector, with names corresponding to markers, and values corresponding to the chromosome on which the named marker is located.
#' @examples
#' map <- qtl::sim.map()
#' pedigree <- f2Pedigree(1000)
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
#' mappedCross <- mpcrossMapped(cross = cross, map = map)
#' chromosomeAssignment <- getChromosomes(mappedCross, markers(mappedCross))
#' chromosomeAssignment
#' @export
getChromosomes <- function(mpcrossMapped, markers)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(missing(markers) || !is.character(markers) || any(is.na(markers)))
	{
		stop("Input markers must be a character vector, without missing values")
	}
	chromosomeIndices <- unlist(sapply(markers, function(marker)
	{
		which(unlist(lapply(mpcrossMapped@map, function(x) !is.na(match(marker, names(x))))))
	}))
	result <- names(mpcrossMapped@map)[chromosomeIndices]
	names(result) <- markers
	return(result)
}
