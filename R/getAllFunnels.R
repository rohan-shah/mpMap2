#' @export
getAllFunnels <- function(cross, standardised = FALSE)
{
	if(!is.logical(standardised) || length(standardised) != 1)
	{
		stop("Input standardised must be TRUE or FALSE")
	}
	if(class(cross) == "geneticData")
	{
		return(.Call("getAllFunnels", cross, PACKAGE="mpMap2"))
	}
	else if(inherits(cross, "mpcross"))
	{
		if(length(cross@geneticData) == 1)
		{
			return(.Call("getAllFunnels", cross@geneticData[[1]], standardised, PACKAGE="mpMap2"))
		}
		else
		{
			return(lapply(cross@geneticData, function(x) .Call("getAllFunnels", x, standardised, PACKAGE="mpMap2")))
		}
	}
	else
	{
		stop("Input must be of class geneticData or mpcross")
	}
}
