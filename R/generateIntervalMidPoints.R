#' @export
generateIntervalMidPoints <- function(object)
{
	result <- lapply(as.list(names(object@map)), function(chrName)
		{
			uniquePositions <- unique(object@map[[chrName]])
			midPoints <- (head(uniquePositions, -1) + tail(uniquePositions, -1))/2
			names(midPoints) <- paste0("Chr", chrName, "Interval", 1:length(positions))
			positions
		})
	names(result) <- names(object@map)
	result
}
