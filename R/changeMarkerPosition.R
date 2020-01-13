#' @title Change the position of a single marker
#' @description Change the position of a single marker
#' @details For an object of class \code{mpcrossMapped}, change the position of a single marker
#' @param mpcrossMapped The object of class \code{mpcrossMapped}, containing a marker to be modified
#' @param marker The name of the marker to change
#' @param newChromosome The new chromosome for the specified marker
#' @param newPosition The new position for the specified marker in cM, on the new chromosome
#' @return A copy of the input object, with the specified marker shifted to the new position and chromosome. 
#' @export
changeMarkerPosition <- function(mpcrossMapped, marker, newChromosome, newPosition)
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(!is.character(marker) || length(marker) != 1 || is.na(marker))
	{
		stop("Input marker must be a marker name")
	}
	if(!is.character(newChromosome) || length(newChromosome) != 1 || is.na(newChromosome) || !(newChromosome %in% names(mpcrossMapped@map)))
	{
		stop("Input newChromosome must be a chromosome name")
	}
	if(!is.numeric(newPosition) || length(newPosition) != 1 || is.na(newPosition) || newPosition < 0)
	{
		stop("Input newPosition must be a marker position")
	}
	originalMarkerOrder <- unlist(lapply(mpcrossMapped@map, names))
	originalPosition <- match(marker, originalMarkerOrder)
	if(is.na(originalPosition)) stop("Could not find marker")

	originalChromosome <- getChromosomes(mpcrossMapped, markers = marker)
	#Remove marker from map
	newMap <- mpcrossMapped@map
	newMap[[originalChromosome]] <- newMap[[originalChromosome]][-match(marker, names(newMap[[originalChromosome]]))]

	#Add marker back in, in new position.
	newMap[[newChromosome]] <- c(newPosition, newMap[[newChromosome]])
	names(newMap[[newChromosome]])[1] <- marker
	newMap[[newChromosome]] <- sort(newMap[[newChromosome]])

	newMarkerOrder <- unlist(lapply(newMap, names))
	newPosition <- match(marker, newMarkerOrder)

	if(newPosition == originalPosition)
	{
		mpcross <- mpcrossMapped
	}
	else
	{
		mpcross <- subset(mpcrossMapped, markers = newMarkerOrder)
	}
	return(new("mpcrossMapped", mpcross, map = newMap))
}
