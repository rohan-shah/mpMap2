setClass("removeHets", contains="NULL")
#' @export
removeHets <- function()
{
	return(new("removeHets"))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "removeHets"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Removing hets will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	for(i in 1:length(e1@geneticData))
	{
		newResults <- .Call("removeHets", e1@geneticData[[i]]@founders, e1@geneticData[[i]]@finals, e1@geneticData[[i]]@hetData, PACKAGE="mpMap2")
		e1@geneticData[[i]]@finals <- newResults$finals
		names(newResults$hetData) <- names(e1@geneticData[[i]]@hetData)
		e1@geneticData[[i]]@hetData <- new("hetData", newResults$hetData)
	}
	return(e1)
})
