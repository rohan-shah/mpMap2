#' @export
addExtraMarker <- function(mpcrossMapped, newMarker, useOnlyExtraImputationPoints = TRUE, intervalMarkerRadius = 50, maxOffset)
{
	if(!require(mpMapInteractive2, quietly = TRUE))
	{
		stop("Function addExtraMarker requires package mpMapInteractiv2")
	}
	if(!inherits(newMarker, "mpcross"))
	{
		stop("Input newMarker must be an object of class mpcross")
	}
	if(!inherits(mpcrossMapped, "mpcrossMapped"))
	{
		stop("Input mpcrossMapped must be an object of class mpcrossMapped")
	}
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Input mpcrossMapped must contain a single data set")
	}
	if(is.null(mpcrossMapped@geneticData[[1]]@imputed))
	{
		stop("Imputation data must be available")
	}
	newMarker <- estimateRF(newMarker, gbLimit = mpcrossMapped@rf@gbLimit, recombValues = mpcrossMapped@rf@theta@levels)
	combined <- combineKeepRF(mpcrossMapped, newMarker, verbose = TRUE, gbLimit = mpcrossMapped@rf@gbLimit, callEstimateRF = TRUE)
	marginalNewMarker <- table(finals(newMarker))
	nObservations <- sum(marginalNewMarker)
	if(useOnlyExtraImputationPoints)
	{
		#Get out the imputation map
		imputationMap <- mpcrossMapped@geneticData[[1]]@imputed@map
		#Flatten the imputation map
		flattenedImputationMapNames <- unlist(lapply(imputationMap, names))
		flattenedImputationMapPositions <- unlist(imputationMap)
		names(flattenedImputationMapPositions) <- flattenedImputationMapNames
		#Get out imputation results only for the grid points
		imputationGridResults <- mpcrossMapped@geneticData[[1]]@imputed@data[,extraImputationPoints(mpcrossMapped)]
		#Compute chi squared statistics
		chiSquared <- apply(imputationGridResults, 2, function(x)
		{
			x[x > 8] <- NA
			subsettedX <- x[!is.na(finals(newMarker))]
			observed <- table(subsettedX, na.omit(finals(newMarker)))
			marginalImputed <- table(subsettedX)
			expected <- outer(marginalImputed, marginalNewMarker) / nObservations
			return(sum((observed - expected)^2 / expected))
		})
		#Name of the best location
		bestLocation <- names(which.max(chiSquared))
		
		chromosomeAssignments <- rep(names(imputationMap), times = unlist(lapply(imputationMap, length)))
		
		#Chromosome for the best location. 
		bestChromosome <- chromosomeAssignments[match(bestLocation, flattenedImputationMapNames)]
		bestPosition <- flattenedImputationMapPositions[bestLocation]
		#Find the index of the marker just to the left
		relevantChromosomeMap <- mpcrossMapped@map[[bestChromosome]]
		markerIndex <- tail(which(relevantChromosomeMap < bestPosition), 1)
		#Get out range. 
		markerRange <- max(1, markerIndex - intervalMarkerRadius):min(length(relevantChromosomeMap), markerIndex + intervalMarkerRadius)

		relevantSubset <- subset(combined, markers = c(names(relevantChromosomeMap)[markerRange], markers(newMarker)))
		grouped <- formGroups(relevantSubset, groups = 1, clusterBy = "theta")
		reorderedGrouped <- mpMapInteractive2(grouped)
		#Get out the reordered markers
		reorderedMarkers <- markers(reorderedGrouped$object)
		permutation <- sapply(reorderedMarkers, function(x) match(x, markers(relevantSubset)))
		if(cor(permutation, 1:length(reorderedMarkers)) < 0) reorderedMarkers <- rev(reorderedMarkers)

		firstMarkerRange <- match(head(markers(grouped), 1), markers(mpcrossMapped))
		lastMarkerRange <- match(tail(markers(grouped), 2)[1], markers(mpcrossMapped))

		#Construct a new overall ordering of *all* the markers (not just the ones on this chromosome)
		newMarkerOrder <- c()
		if(firstMarkerRange != 1) newMarkerOrder <- markers(mpcrossMapped)[1:(firstMarkerRange-1)]

		newMarkerOrder <- c(newMarkerOrder, reorderedMarkers)
		if(lastMarkerRange != nMarkers(mpcrossMapped)) newMarkerOrder <- c(newMarkerOrder, markers(mpcrossMapped)[(lastMarkerRange+1):nMarkers(mpcrossMapped)])

		#Put original object (Plus extra marker) in the new order. 
		objectInNewOrder <- subset(combined, markers = newMarkerOrder)
		newChromosomeMarkers <- newMarkerOrder[(newMarkerOrder %in% names(mpcrossMapped@map[[bestChromosome]])) | newMarkerOrder == markers(newMarker)]
		#Get out just the chromosome that's been changed. We only need to recompute the map for these markers.
		changedChromosomeInNewOrder <- subset(objectInNewOrder, markers = newChromosomeMarkers)

		#Estimate RF fractions for the changed chromosome
		changedChromosomeInNewOrder <- formGroups(changedChromosomeInNewOrder, groups = 1, clusterBy = "theta")
		#Re-estimate map for the changed chromosome
		newMap <- estimateMap(changedChromosomeInNewOrder, maxOffset = maxOffset, mapFunction = haldane, verbose = TRUE)
		#Copy the old map, and update the chromosome that changed. 
		finalMap <- mpcrossMapped@map
		finalMap[[bestChromosome]] <- newMap[[1]]
		#Put the new map into the new object
		objectInNewOrder <- new("mpcrossMapped", objectInNewOrder, map = finalMap, rf = objectInNewOrder@rf)
		#Check that we have a valid object still.
		validObject(objectInNewOrder, complete = TRUE)
		return(list(statistics = chiSquared, object = objectInNewOrder))
	}
	else
	{
		stop("This code path is not implemented yet")
	}
}
