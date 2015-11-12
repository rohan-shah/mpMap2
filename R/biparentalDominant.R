setClass("biparentalDominant", contains="NULL")
#' @export
biparentalDominant <- function()
{
	return(new("biparentalDominant"))
}
setMethod(f = "+", signature = c("geneticData", "biparentalDominant"), definition = function(e1, e2)
{
	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		hetDataThisMarker <- copied@hetData[[x]]
		isHetRow <- hetDataThisMarker[,1] != hetDataThisMarker[,2]
		hetCode <- unique(hetDataThisMarker[isHetRow, 3])
		homozygoteCodes <- hetDataThisMarker[!isHetRow, 3]
		dominant <- sample(homozygoteCodes, 1)

		copied@hetData[[x]][isHetRow,3] <<- dominant
		copied@finals[copied@finals[,x] == hetCode,x] <<- dominant
	})
	return(copied)
})
setMethod(f = "+", signature = c("mpcross", "biparentalDominant"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Assigning random dominant marker patterns will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	if(length(e1@geneticData) > 1)
	{
		stop("Attempting to change an object containing multiple data sets. Please change each dataset individually")
	}
	if(nFounders(e1) != 2)
	{
		stop("biparentalDomimant can only be applied to a biparental design")
	}
	e1@geneticData[[1]] <- e1@geneticData[[1]]+e2
	return(e1)
})
