expand <- function(mpcross, newMarkers)
{
	inheritsNewMpcrossArgument(mpcross)
	if(class(mpcross) != "mpcross")
	{
		warning(paste0("Converting object of class ", class(mpcross), " to class mpcross, data may be lost"))
		mpcross <- as(mpcross, "mpcross")
	}

	oldMarkers <- markers(mpcross)
	nOldMarkers <- length(oldMarkers)
	nNewMarkers <- length(newMarkers)
	if(!all(oldMarkers %in% newMarkers))
	{
		stop("New marker set must contain old marker set")
	}

	newGeneticData <- lapply(mpcross@geneticData, function(x)
	{
		newFounders <- matrix(data = NA, nrow = nFounders(x), ncol = nNewMarkers)
		newFinals <- matrix(data = NA, nrow = nLines(x), ncol = nNewMarkers)
		colnames(newFounders) <- colnames(newFinals) <- newMarkers
		rownames(newFounders) <- rownames(x@founders)
		rownames(newFinals) <- rownames(x@finals)

		newFounders[,oldMarkers] <- x@founders
		newFinals[,oldMarkers] <- x@finals

		newHetData <- replicate(nNewMarkers, matrix(0L, 0, 3), simplify=NULL)
		names(newHetData) <- newMarkers
		newHetData[oldMarkers] <- x@hetData

		return(new("geneticData", finals = newFinals, founders = newFounders, pedigree = x@pedigree, hetData = new("hetData", newHetData)))
	})
	return(new("mpcross", geneticData = newGeneticData))
}