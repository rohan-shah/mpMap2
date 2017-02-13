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
