fullHetData <- function(map, nFounders)
{
	nonNegativeIntegerArgument(nFounders)
	isMapArgument(map)

	repeatedValue <- cbind(as.matrix(expand.grid(1:nFounders, 1:nFounders)), 1:(nFounders*nFounders))
	dimnames(repeatedValue) <- NULL
	allMarkers <- unlist(lapply(map, names))
	retVal <- replicate(length(allMarkers), repeatedValue, simplify=FALSE)
	names(retVal) <- allMarkers
	return(new("hetData", retVal))
}