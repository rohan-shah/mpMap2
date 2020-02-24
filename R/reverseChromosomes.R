#' @title Reverse the order of the specified chromosomes
#' @description Create a new object, with the specified chromosomes reversed
#' @details Create a new object, with the specified chromosomes reversed
#' @param mpcrossMapped The initial object, for which we want to reverse some of the chromosomes
#' @param chromosomes The names of the chromosomes to reverse
#' @return An object of class \code{mpcrossMapped}, with certain chromosomes reversed.
#' @examples
#' map <- qtl::sim.map()
#' pedigree <- f2Pedigree(1000)
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
#' mappedCross <- mpcrossMapped(cross = cross, map = map)
#' reversedX <- reverseChromosomes(mappedCross, "X")
#' reversedX@map[["X"]]
#' mappedCross@map[["X"]]
#' @export
reverseChromosomes <- function(mpcrossMapped, chromosomes)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(!all(chromosomes %in% names(mpcrossMapped@map)))
	{
		stop("Invalid values in input chromosomes")
	}
	newMap <- mpcrossMapped@map
	for(chromosome in chromosomes)
	{
		newMap[[chromosome]] <- rev(newMap[[chromosome]])
		newMap[[chromosome]] <- max(newMap[[chromosome]]) - newMap[[chromosome]]
	}
	if(any(unlist(lapply(mpcrossMapped@geneticData, function(x) !is.null(x@probabilities)))))
	{
		warning("Discarding probability data")
	}
	if(any(unlist(lapply(mpcrossMapped@geneticData, function(x) !is.null(x@imputed)))))
	{
		warning("Discarding imputation data")
	}
	newOrder <- unlist(lapply(newMap, names))
	newRF <- NULL
	if(!is.null(mpcrossMapped@rf))
	{
		newRF <- subset(mpcrossMapped@rf, markers = newOrder)
	}
	return(new("mpcrossMapped", subset(mpcrossMapped, markers = newOrder), map = newMap, rf = newRF))
}
