expand <- function(mpcross, newMarkers)
{
	inheritsNewMpcrossArgument(mpcross)

	oldMarkers <- markers(mpcross)
	nOldMarkers <- length(oldMarkers)
	nNewMarkers <- length(newMarkers)
	if(!all(oldMarkers %in% newMarkers))
	{
		stop("New marker set must contain old marker set")
	}
	#If the new markers and the old markers are the same, just return the object
	if(all(newMarkers %in% oldMarkers))
	{
		return(subset(mpcross, markers = newMarkers))
	}

	if(class(mpcross) != "mpcross")
	{
		#Warn that RF data will be lost
		if(inherits(mpcross, "mpcrossRF"))
		{
			warning(paste0("Converting object of class mpcrossRF to class mpcross, recombination data will be lost"))
		}
		if(inherits(mpcross, "mpcrossLG"))
		{
			if(!is.null(mpcross@rf))
			{
				warning(paste0("Converting object of class mpcrossLG to class mpcross, recombination data will be lost"))
			}
			warning(paste0("Converting object of class mpcrossLG to class mpcross, linkage group data will be lost"))
		}
		mpcross <- as(mpcross, "mpcross")
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