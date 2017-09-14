#' @export
generateGridPositions <- function(spacing)
{
	retFunction <- function(object)
	{
		result <- lapply(as.list(names(object@map)), function(chrName)
			{
				x <- object@map[[chrName]]
				range <- range(x)
				positions <- seq(range[1], range[2], by = spacing)
				names(positions) <- paste0("Chr", chrName, "Loc", 1:length(positions))
				positions
			})
		names(result) <- names(object@map)
		result
	}
	return(retFunction)
}
