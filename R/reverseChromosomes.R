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
	return(new("mpcrossMapped", subset(mpcrossMapped, markers = newOrder), map = newMap))
}
