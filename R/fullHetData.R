fullHetData <- function(map, nFounders)
{
	nonNegativeIntegerArgument(nFounders)
	isMapArgument(map)

	repeatedValue <- cbind(as.matrix(expand.grid(1:nFounders, 1:nFounders)), 0L)
	repeatedValue[seq(1, nFounders*nFounders, by = nFounders+1),3] <- 1:nFounders

	#Put in allele values
	for(row in 1:nrow(repeatedValue))
	{
		if(repeatedValue[row, 3] == 0)
		{
			if(repeatedValue[row, 2] > repeatedValue[row, 1])
			{
				repeatedValue[row, 3] <- repeatedValue[nFounders * (repeatedValue[row,1]-1) + repeatedValue[row,2],3]
			}
			else repeatedValue[row, 3] <- max(repeatedValue[,3])+1L
		}
	}
	dimnames(repeatedValue) <- NULL
	allMarkers <- unlist(lapply(map, names))
	retVal <- replicate(length(allMarkers), repeatedValue, simplify=FALSE)
	names(retVal) <- allMarkers
	return(new("hetData", retVal))
}