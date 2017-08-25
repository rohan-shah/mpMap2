splitVector <- function(vector, splitValue)
{
	index <- match(splitValue, vector)
	if(index == length(vector)) return(list(before = vector[1:length(vector)], after = c()))
	return(list(before = vector[1:index], after = vector[(index+1):length(vector)]))
}
#' @export
addExtraMarkers <- function(mpcrossMapped, newMarkers, useOnlyExtraImputationPoints = TRUE, reorderRadius = 103, maxOffset = 50, knownChromosome, imputationArgs = NULL, onlyStatistics = FALSE, orderCrossArgs = list(), attemptMpMap2Interactive = TRUE, verbose = TRUE, reorder = TRUE)
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
	if(reorderRadius < 2*maxOffset + 3)
	{
		stop("Input reorderRadius must be at least 2*maxOffset + 3")
	}
	newMarkers <- estimateRF(newMarkers, gbLimit = mpcrossMapped@rf@gbLimit, recombValues = mpcrossMapped@rf@theta@levels)
	founders <- nFounders(mpcrossMapped)
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
		finalSubset <- finals(newMarkers)[,1]
		finalSubsetIsNA <- is.na(finalSubset)
		finalSubsetOmitNA <- na.omit(finalSubset)
		#Compute chi squared statistics
		chiSquared <- apply(imputationGridResults, 2, function(x)
		{
			x[x > founders] <- NA
			subsettedX <- x[!finalSubsetIsNA]
			observed <- table(subsettedX, finalSubsetOmitNA)
			marginalImputed <- table(subsettedX)
			marginalNewMarker <- apply(observed, 2, sum)
			expected <- outer(marginalImputed, marginalNewMarker) / sum(marginalImputed)
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
		markerIndex <- tail(which(relevantChromosomeMap <= bestPosition), 1)
		#Get out range. 
		markerRange <- max(1, markerIndex - reorderRadius):min(length(relevantChromosomeMap), markerIndex + reorderRadius)
		cat("Reordering markers [", min(markerRange), ":", max(markerRange), "] of the chromosome, total nmuber of markers for this chromosome was ", length(relevantChromosomeMap), "\n", sep = "")

		#Combine and keep recombination fractions
		combined <- combineKeepRF(mpcrossMapped, newMarkers, verbose = verbose, gbLimit = mpcrossMapped@rf@gbLimit, callEstimateRF = match(names(mpcrossMapped@map[[bestChromosome]]), markers(mpcrossMapped)))
		#Put the new subset of markers in the right place.
		splitResults <- splitVector(names(relevantChromosomeMap)[markerRange], names(relevantChromosomeMap)[markerIndex])
		relevantSubset <- subset(combined, markers = c(splitResults$before, markers(newMarkers), splitResults$after))
		grouped <- formGroups(relevantSubset, groups = 1, clusterBy = "theta")
		if(reorder)
		{
			if(hasMpMapInteractive2)
			{
				reorderedGrouped <- mpMapInteractive2(grouped)$object
			}
			else
			{
				reorderedGrouped <- do.call(orderCross, c(list(grouped), orderCrossArgs))
			}
		}
		else reorderedGrouped <- grouped
		#Get out the reordered markers
		reorderedMarkers <- markers(reorderedGrouped)

		permutation <- sapply(reorderedMarkers, function(x) match(x, markers(relevantSubset)))
		if(cor(permutation, 1:length(reorderedMarkers)) < 0) reorderedMarkers <- rev(reorderedMarkers)

		firstMarkerRange <- match(names(relevantChromosomeMap)[min(markerRange)], markers(mpcrossMapped))
		lastMarkerRange <- match(names(relevantChromosomeMap)[max(markerRange)], markers(mpcrossMapped))

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
		if(all(reorderedMarkers == markers(grouped)))
		{
			#If we didn't reorder anything interactively, then do the update quicker. 
			newMap <- estimateMap(reorderedGrouped, maxOffset = maxOffset, mapFunction = haldane, verbose = verbose)
			names(newMap) <- bestChromosome
			#Copy the old map, and update *part* of the chromosome that changed. 
			finalMap <- mpcrossMapped@map
			#Only update the bit in the middle, with a radius of maxOffset. 
			smallerMarkerRange <- max(1, markerIndex - maxOffset):min(length(relevantChromosomeMap), markerIndex + maxOffset)
			#Start of the added chunk of markers, within the new map chunk
			startNewRangeWithin <- match(head(markers(newMarkers), 1), reorderedMarkers)
			#End of the added chunk of markers, within the new map chunk
			endNewRangeWithin <- match(tail(markers(newMarkers), 1), reorderedMarkers)
			#Range of positions within the re-estimated map, to use
			withinReestimatedMapRange <- max(1, startNewRangeWithin - maxOffset):min(length(newMap[[1]]), endNewRangeWithin + maxOffset)
			#First bit stays the same. Unless the reestimated bit covers the first marker, in which case the new map *starts* with the re-estimated bit. 
			if(min(withinReestimatedMapRange) == 1) updatedChromosomeMap <- c()
			else updatedChromosomeMap <- finalMap[[bestChromosome]][names(relevantChromosomeMap)[setdiff(1:markerIndex, smallerMarkerRange+1)]]
			#Then the updated bit
			updatedChromosomeMap <- c(updatedChromosomeMap, newMap[[1]][withinReestimatedMapRange] + finalMap[[bestChromosome]][head(names(newMap[[1]][withinReestimatedMapRange]), 1)] - head(newMap[[1]][withinReestimatedMapRange], 1))
			#Then for the last bit the distances between adjacent markers stays the same, but the offset has to change. 
			finalPart <- finalMap[[bestChromosome]][setdiff(names(relevantChromosomeMap)[(markerIndex+1):length(relevantChromosomeMap)], names(relevantChromosomeMap)[smallerMarkerRange])]
			firstUnchanged <- head(names(finalPart), 1)
			finalPart <- finalPart - finalPart[1] + tail(updatedChromosomeMap, 1) + finalMap[[bestChromosome]][firstUnchanged] - finalMap[[bestChromosome]][match(firstUnchanged, names(finalMap[[bestChromosome]])) - 1]
			updatedChromosomeMap <- c(updatedChromosomeMap, finalPart)
			#if(any(order(updatedChromosomeMap) != 1:length(updatedChromosomeMap)))
			#{
			#	stop("Internal error")
			#}
			finalMap[[bestChromosome]] <- updatedChromosomeMap
		}
		else
		{
			newMap <- estimateMap(changedChromosomeInNewOrder, maxOffset = maxOffset, mapFunction = haldane, verbose = verbose)
			names(newMap) <- bestChromosome
			#Copy the old map, and update the chromosome that changed. 
			finalMap <- mpcrossMapped@map
			finalMap[[bestChromosome]] <- newMap[[1]]
		}
		#Put the new map into the new object
		objectInNewOrder <- new("mpcrossMapped", objectInNewOrder, map = finalMap, rf = objectInNewOrder@rf)

		#Update the imputation data, so we don't have to run the entire thing all over again. 
		if(!is.null(imputationArgs))
		{
			previousKey <- mpcrossMapped@geneticData[[1]]@imputed@key
			mappedChangedChromosomeInNewOrder <- new("mpcrossMapped", changedChromosomeInNewOrder, rf = changedChromosomeInNewOrder@rf, map = finalMap[bestChromosome])
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
			#validObject(objectInNewOrder, complete = TRUE, test = TRUE)
		}
		return(list(statistics = chiSquared, object = objectInNewOrder))
	}
	else
	{
		stop("This code path is not implemented yet")
	}
}
