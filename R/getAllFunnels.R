#' @export
getAllFunnels <- function(cross)
{
	if(class(cross) == "geneticData")
	{
		return(.Call("getAllFunnels", cross, PACKAGE="mpMap2"))
	}
	else if(inherits(cross, "mpcross"))
	{
		return(lapply(cross@geneticData, function(x) .Call("getAllFunnels",x, PACKAGE="mpMap2")))
	}
	else
	{
		stop("Input must be of class geneticData or mpcross")
	}
}
