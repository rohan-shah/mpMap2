#' @title Specify an equally spaced grid of genetic positions 
#' @description Specify an equally spaced grid of genetic positions 
#' @details Some functions, such as \code{imputeFounders} and \code{computeGenotypeProbabilities}, take in a set of genetic positions as one of the inputs. This function is an easy way to speify an equally spaced grid of positions. 
#' 
#' Note that the return value is itself a function, which is applied internally by \code{imputeFounders} or \code{computeGenotypeProbabilities} to an object of class \code{mpcrossMapped}. 
#' @param spacing The spacing of the genetic positions, in cM.
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
