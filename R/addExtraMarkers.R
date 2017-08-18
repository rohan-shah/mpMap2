#' @export
addExtraMarkers <- function(mpcrossMapped, newMarkers, useOnlyExtraImputationPoints = TRUE, reorderRadius = 50, maxOffset = 50, knownChromosome, imputationArgs = NULL, onlyStatistics = FALSE, orderCrossArgs = list(), attemptMpMap2Interactive = TRUE, verbose = TRUE)
{
	hasMpMapInteractive2 <- attemptMpMap2Interactive && require(mpMapInteractive2, quietly = TRUE)
	if(!inherits(newMarkers, "mpcross"))
	{
		stop("Input newMarkers must be an object of class mpcross")
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
	if(!missing(knownChromosome) && !(knownChromosome %in% names(mpcrossMapped@map)))
	{
		stop("Input knownChromosome must be a chromosome named in mpcrossMapped@map")
	}
	if(inherits(newMarkers, "mpcrossMapped"))
	{
		warning("Discarding map for additional markers")
		newMarkers <- as(newMarkers, "mpcross")
	}
	newMarkers <- estimateRF(newMarkers, gbLimit = mpcrossMapped@rf@gbLimit, recombValues = mpcrossMapped@rf@theta@levels)
	marginalNewMarker <- table(finals(newMarkers)[,1])
	nObservations <- sum(marginalNewMarker)
	if(useOnlyExtraImputationPoints)
	{
		if(length(extraImputationPoints(mpcrossMapped)) == 0)
		{
			stop("There were no additional imputation points")
		}
		#Get out the imputation map
		imputationMap <- imputationMap(mpcrossMapped)
		#Flatten the imputation map
		flattenedImputationMapNames <- flatImputationMapNames(mpcrossMapped)

		flattenedImputationMapPositions <- unlist(imputationMap)
		names(flattenedImputationMapPositions) <- flattenedImputationMapNames
		#Get out imputation results only for the grid points
		imputationGridResults <- imputationData(mpcrossMapped)[,extraImputationPoints(mpcrossMapped)]
		#Compute chi squared statistics
		chiSquared <- apply(imputationGridResults, 2, function(x)
		{
			x[x > 8] <- NA
			subsettedX <- x[!is.na(finals(newMarkers)[,1])]
			observed <- table(subsettedX, na.omit(finals(newMarkers)[,1]))
			marginalImputed <- table(subsettedX)
			expected <- outer(marginalImputed, marginalNewMarker) / nObservations
			return(sum((observed - expected)^2 / expected))
		})
		if(onlyStatistics)
		{
			return(chiSquared)
		}
		chromosomeAssignments <- rep(names(imputationMap), times = unlist(lapply(imputationMap, length)))
		#Name of the best location
		if(missing(knownChromosome))
		{
			bestLocation <- names(which.max(chiSquared))
		}
		else
		{
			bestLocation <- names(which.max(chiSquared[intersect(names(imputationMap[[knownChromosome]]), extraImputationPoints(mpcrossMapped))]))
		}
		#Chromosome for the best location. 
		bestChromosome <- chromosomeAssignments[match(bestLocation, flattenedImputationMapNames)]

		bestPosition <- flattenedImputationMapPositions[bestLocation]
		#Find the index of the marker just to the left
		relevantChromosomeMap <- mpcrossMapped@map[[bestChromosome]]
		markerIndex <- tail(which(relevantChromosomeMap < bestPosition), 1)
		#Get out range. 
		markerRange <- max(1, markerIndex - reorderRadius):min(length(relevantChromosomeMap), markerIndex + reorderRadius)
		cat("Reordering markers [", min(markerRange), ":", max(markerRange), "] of the chromosome, total nmuber of markers for this chromosome was ", length(relevantChromosomeMap), "\n", sep = "")

		#Combine and keep recombination fractions
		combined <- combineKeepRF(mpcrossMapped, newMarkers, verbose = verbose, gbLimit = mpcrossMapped@rf@gbLimit, callEstimateRF = match(names(mpcrossMapped@map[[bestChromosome]]), markers(mpcrossMapped)))
		#Put the new subset of markers in the right place.
		if(min(markerRange) == 1 && max(markerRange) == length(relevantChromosomeMap))
		{
			relevantSubset <- subset(combined, markers = c(names(relevantChromosomeMap)[1:markerIndex], markers(newMarkers), names(relevantChromosomeMap)[(markerIndex+1):length(relevantChromosomeMap)]))
		}
		else if(min(markerRange) == 1)
		{
			relevantSubset <- subset(combined, markers = c(names(relevantChromosomeMap)[1:markerIndex], markers(newMarkers), names(relevantChromosomeMap)[(markerIndex+1):(markerIndex+reorderRadius)]))
		}
		else if(max(markerRange) == length(relevantChromosomeMap))
		{
			relevantSubset <- subset(combined, markers = c(names(relevantChromosomeMap)[(markerIndex - reorderRadius):markerIndex], markers(newMarkers), names(relevantChromosomeMap)[(markerIndex+1):length(relevantChromosomeMap)]))
		}
		else
		{
			relevantSubset <- subset(combined, markers = c(names(relevantChromosomeMap)[(markerIndex - reorderRadius):markerIndex], markers(newMarkers), names(relevantChromosomeMap)[(markerIndex+1):(markerIndex+reorderRadius)]))
		}
		grouped <- formGroups(relevantSubset, groups = 1, clusterBy = "theta")
		if(hasMpMapInteractive2)
		{
			reorderedGrouped <- mpMapInteractive2(grouped)
			#Get out the reordered markers
			reorderedMarkers <- markers(reorderedGrouped$object)
		}
		else
		{
			reorderedGrouped <- do.call(orderCross, c(list(grouped), orderCrossArgs))
			#Get out the reordered markers
			reorderedMarkers <- markers(reorderedGrouped)
		}
		permutation <- sapply(reorderedMarkers, function(x) match(x, markers(relevantSubset)))
		if(cor(permutation, 1:length(reorderedMarkers)) < 0) reorderedMarkers <- rev(reorderedMarkers)

		firstMarkerRange <- match(names(relevantChromosomeSubset)[min(markerRange)], markers(mpcrossMapped))
		lastMarkerRange <- match(names(relevantChromosomeSubset)[max(markerRange)], markers(mpcrossMapped))

		#Construct a new overall ordering of *all* the markers (not just the ones on this chromosome)
		newMarkerOrder <- c()
		if(firstMarkerRange != 1) newMarkerOrder <- markers(mpcrossMapped)[1:(firstMarkerRange-1)]

		newMarkerOrder <- c(newMarkerOrder, reorderedMarkers)
		if(lastMarkerRange != nMarkers(mpcrossMapped)) newMarkerOrder <- c(newMarkerOrder, markers(mpcrossMapped)[(lastMarkerRange+1):nMarkers(mpcrossMapped)])

		#Put original object (Plus extra marker) in the new order. 
		objectInNewOrder <- subset(combined, markers = newMarkerOrder)
		newChromosomeMarkers <- newMarkerOrder[(newMarkerOrder %in% names(mpcrossMapped@map[[bestChromosome]])) | (newMarkerOrder %in% markers(newMarkers))]
		#Get out just the chromosome that's been changed. We only need to recompute the map for these markers.
		changedChromosomeInNewOrder <- subset(objectInNewOrder, markers = newChromosomeMarkers)

		#Estimate RF fractions for the changed chromosome
		changedChromosomeInNewOrder <- formGroups(changedChromosomeInNewOrder, groups = 1, clusterBy = "theta")
		#Re-estimate map for the changed chromosome
		newMap <- estimateMap(changedChromosomeInNewOrder, maxOffset = maxOffset, mapFunction = haldane, verbose = verbose)
		names(newMap) <- bestChromosome
		#Copy the old map, and update the chromosome that changed. 
		finalMap <- mpcrossMapped@map
		finalMap[[bestChromosome]] <- newMap[[1]]
		#Put the new map into the new object
		objectInNewOrder <- new("mpcrossMapped", objectInNewOrder, map = finalMap, rf = objectInNewOrder@rf)

		#Update the imputation data, so we don't have to run the entire thing all over again. 
		if(!is.null(imputationArgs))
		{
			previousKey <- mpcrossMapped@geneticData[[1]]@imputed@key
			mappedChangedChromosomeInNewOrder <- new("mpcrossMapped", changedChromosomeInNewOrder, rf = changedChromosomeInNewOrder@rf, map = newMap)
			changedChromosomeInNewOrder <- do.call(imputeFounders, c(list(mappedChangedChromosomeInNewOrder), imputationArgs))
			if(!identical(previousKey, changedChromosomeInNewOrder@geneticData[[1]]@imputed@key))
			{
				stop("Generated imputation data was incompatible with previous imputation data; could not update old data")
			}
			#Form the new imputation map
			newImputationMap <- imputationMap
			newImputationMap[[bestChromosome]] <- imputationMap(changedChromosomeInNewOrder)[[1]]
			
			#Construct the new imputation data matrix
			newImputationData <- matrix(-1L, nrow = nLines(objectInNewOrder), ncol = sum(unlist(lapply(newImputationMap, length))))
			rownames(newImputationData) <- lineNames(objectInNewOrder)
			colnames(newImputationData) <- unlist(lapply(newImputationMap, names))

			#Put in the data. First the original data
			oldMapWithoutBestChromosome <- imputationMap
			oldMapWithoutBestChromosome[[bestChromosome]] <- NULL
			oldMapWithoutBestChromosomeMarkers <- unlist(lapply(oldMapWithoutBestChromosome, names))
			newImputationData[,oldMapWithoutBestChromosomeMarkers] <- imputationData(mpcrossMapped)[,oldMapWithoutBestChromosomeMarkers]
			#And then the recomputed data for the chromosome that changed. 
			newImputationData[,flatImputationMapNames(changedChromosomeInNewOrder)] <- imputationData(changedChromosomeInNewOrder)

			#Form the new imputation data object
			newImputationObject <- new("imputed", data = newImputationData, key = previousKey, map = newImputationMap, errors = NULL)
			#Insert into overall object
			objectInNewOrder@geneticData[[1]]@imputed <- newImputationObject
			validObject(objectInNewOrder, complete = TRUE, test = TRUE)
		}
		return(list(statistics = chiSquared, object = objectInNewOrder))
	}
	else
	{
		stop("This code path is not implemented yet")
	}
}
