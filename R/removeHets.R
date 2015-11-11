setClass("removeHets", contains="NULL")
#' @export
removeHets <- function()
{
	return(new("removeHets"))
}
setMethod(f = "+", signature = c("mpcross", "removeHets"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Removing hets will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	lapply(1:length(e1@geneticData), 
		function(index)
		{
			founders <- nFounders(e1@geneticData[[index]])
			e1@geneticData[[index]]@finals[e1@geneticData[[index]]@finals > founders] <<- NA
		})
	return(e1)
})
