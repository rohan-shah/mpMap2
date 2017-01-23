#' @export
testDistortion <- function(object)
{
	if(inherits(object, "mpcross") && length(object@geneticData) == 1)
	{
		object <- object@geneticData[[1]]
	}
	if(!isS4(object) || !inherits(object, "geneticData"))
	{
		stop("Input object must be an S4 object of class geneticData")
	}
	result <- .Call("testDistortion", object, PACKAGE="mpMap2")
	pValues <- 1 - pchisq(result$testStatistics, result$classes - 1)
	result <- data.frame(pValue = pValues, statistic = result$testStatistics, L1 = result$L1, L2 = result$L2)
	rownames(result) <- markers(object)
	return(result)
}
