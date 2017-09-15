#' @export
generateIntervalMidPoints <- function(object)
{
	if(!is(object, "mpcrossMapped"))
	{
		stop("Input must be an object of class mpcrossMapped")
	}
	result <- lapply(as.list(names(object@map)), function(chrName)
		{
			uniquePositions <- unique(object@map[[chrName]])
			midPoints <- (head(uniquePositions, -1) + tail(uniquePositions, -1))/2
			names(midPoints) <- paste0("Chr", chrName, "Interval", 1:length(midPoints))
			midPoints
		})
	names(result) <- names(object@map)
	result
}
