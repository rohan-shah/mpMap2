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
