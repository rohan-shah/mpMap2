#' @include map-class.R
splitVector <- function(vector, splitValue)
{
	index <- match(splitValue, vector)
	if(index == length(vector)) return(list(before = vector[1:length(vector)], after = c()))
	return(list(before = vector[1:index], after = vector[(index+1):length(vector)]))
}
.addExtraMarkersStatistics <- setClass("addExtraMarkersStatistics", slots=list(data = "numeric", map = "map"))
#' @title Plot chi-squared statistics for independence
#' @description Plot the chi-squared statistics for independence, used to map a new marker to an existing genetic map
#' @details This function plots a trace of the chi-squared test-statistics used to map a new genetic marker to an existing genetic map. This can be useful to, for example, see if a single polymorphism is present at multiple points on the genome. 
#' @param x Object of class \code{addExtraMarkersStatistics} containing test-statistic values.
#' @param y Unused
#' @param ... Unused
#' @return A ggplot object suitable for display. 
setMethod(f = "plot", signature = "addExtraMarkersStatistics", definition = function(x, ...)
{
	offsetVector <- c(0, head(unlist(lapply(x@map, max)), -1))
	cumsumOffsetVector <- cumsum(offsetVector)
	offsetVector <- rep(cumsumOffsetVector, times = unlist(lapply(x@map, length)))
	data <- data.frame(position = offsetVector + unlist(x@map), statistic = x@data)
	verticalLineData <- data.frame(position = tail(cumsumOffsetVector, -1))
	return(ggplot(data, mapping = aes_string(x = 'position', y = 'statistic')) + geom_line(color = "red") + xlab("Distance (cM)") + ylab("Test statistic") + geom_vline(data = verticalLineData, mapping = aes_string(xintercept = 'position')))
})
#' @title Add extra markers
#' @description Add extra markers to a map, using a QTL-mapping style approach. 
#' @param mpcrossMapped An object of class \code{mpcrossMapped} (dataset with a map), which must include imputed IBD genotypes and recombination fraction data. 
#' @param newMarkers An object of class \code{mpcross} containing the new markers to add. 
#' @param useOnlyExtraImputationPoints Should we only attempt to add the new marker at points at which imputation data has been generated, which are \emph{not} markers? Currently this must be \code{TRUE}. In future \code{FALSE} may be allowed. 
#' @param reorderRadius The width of the region (in terms of number of markers) in which to attempt to reorder, after the extra markers are added. 
#' @param maxOffset The maxOffset parameter for the call to \code{\link{estimateMap}}, which is used to re-estimate the map (locally), after the additional markers are added. 
#' @param knownChromosome The name of a chromosome, if the extra markers are known to go on a specific chromosome
#' @param imputationArgs A list containing additional arguments to \code{\link{imputeFounders}}. 
#' @param onlyStatistics If this argument is \code{TRUE}, then only the chi-squared test statistic values are computed, and the markers are not added.
#' @param orderCrossArgs A list containing additional arguments to \code{\link{orderCross}}. 
#' @param verbose Should extra logging output be generated?
#' @param reorder Should local reordering be performed after the extra markers are added?
#' @return If \code{onlyStatistics} was set to \code{TRUE}, an object of class \code{addExtraMarkersStatistics} containing the test statistic values. If \code{onlyStatistics} was set to \code{FALSE}, a list containing the test statistic values in entry \code{statistics} and in entry \code{object}, a new object containing the input object with the new markers added.
#' 
#' @details This function uses a QTL-mapping style approach to add extra markers to an existing map. The code uses the imputation data at a collection of points, and the marker alleles for the \emph{first} marker of the extra markers. If the imputed genotype at a point is \emph{independent} from the genotype at the new marker, then the new marker probably should \emph{not} be mapped to that point. If the imputed genotype at a point and the marker allele are \emph{strongly dependent}, then the new marker \emph{should} probably be mapped to that point. Dependence and independence are measured using a chi-squared test stastistic for independence. \emph{All the extra markers} are then mapped to the position where the test statistic is largest. It is recommended that only single markers be added at a time, unless you are extremely confident that all the extra markers should be located at the same position.
#'
#' Currently the set of points at which the new markers are considered for addition is the set of points at which imputation data is available, \emph{which are not markers}. The intention is that this set of points should be an equally spaced grid of points; this reduced the number of tests that are performed, as generally there are far fewer points in the grid, than there are markers. After the new marker is added, local reordering will need to be performed anyway, making any loss in accuracy by using the grid of points largely irrelevant. In future it may be possible to use the set of all marker positions as the set of points at which tests are performed, by setting \code{useOnlyExtraImputationPoints} to \code{FALSE}. 
#' 
#' Once the extra markers have added, local reordering is optionally performed, depending on argument \code{reordering}. The radius of the region on which reordering is performed, in terms of the number of markers, is \code{reorderRadius}. 
#' 
#' Once the optional reordering step has been performed, the map is recomputed locally, to account for the addition of the extra markers. The argument \code{maxOffset} is passed through to \code{estimateMap}. Finally, the imputation data will be recomputed if \code{imputationArgs} is not \code{NULL}; in that case, \code{imputationArgs} should contain a list of arguments to \code{imputeFounders}. It is recommended that the imputation data be recomputed if further markers are to be added. 
#' 
#' In some cases the user will want to apply a threshold to the maximum value of the test statistics, and only add the marker if the test statistics exceed the threshold. In this case the function should be called twice. For the first call, \code{onlyStatistics} should be set to \code{FALSE}. If the resulting test statistics exceed the threshold, then \code{addExtrMarkers} should be called again with \code{onlyStatistics} set to \code{TRUE}. 
#'
#' @examples
#' data(simulatedFourParentData)
#' #Create object that includes the correct map
#' mapped <- new("mpcrossMapped", simulatedFourParentData, map = simulatedFourParentMap)
#' #Remove marker number 50. Normally the map is discarded, but we specify to keep it. 
#' removedMiddle <- subset(mapped, markers = (1:101)[-50], keepMap = TRUE)
#' #Compute imputation data, at all the markers, and an equally spaced grid of points
#' removedMiddle <- imputeFounders(removedMiddle, errorProb = 0.02, 
#' 	extraPositions = generateGridPositions(1))
#' #Estimate recombination fractions
#' removedMiddle <- estimateRF(removedMiddle)
#' #Get out the extra marker to add
#' extraMarker <- subset(simulatedFourParentData, markers = 50)
#' #Add the extra marker, without doing any local reordering. After the marker is added, 
#' #	recompute the imputation data, using the same arguments as previously. 
#' \dontrun{withExtra <- addExtraMarkers(mpcrossMapped = removedMiddle, newMarkers = extraMarker, 
#' 	reorder = FALSE, imputationArgs = list(errorProb = 0.02, 
#' 	extraPositions = generateGridPositions(1)))$object}
#' @export
addExtraMarkers <- function(mpcrossMapped, newMarkers, useOnlyExtraImputationPoints = TRUE, reorderRadius = 103, maxOffset = 50, knownChromosome, imputationArgs = NULL, onlyStatistics = FALSE, orderCrossArgs = list(), 
#attemptMpMap2Interactive = TRUE, 
			    verbose = TRUE, reorder = TRUE)
{
#	hasMpMapInteractive2 <- attemptMpMap2Interactive && requireNamespace("mpMapInteractive2", quietly = TRUE)
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
	if(is.null(mpcrossMapped@rf))
	{
		stop("Input mpcrossMapped must include recombination fraction data")
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
	if(!onlyStatistics)
	{
		newMarkers <- estimateRF(newMarkers, gbLimit = mpcrossMapped@rf@gbLimit, recombValues = mpcrossMapped@rf@theta@levels)
	}
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
		if(length(unique(finalSubset)) == 1)
		{
			stop("Cannot map a marker which is monomorphic")
		}
		finalSubsetIsNA <- is.na(finalSubset)
		finalSubsetOmitNA <- na.omit(finalSubset)
		#Compute chi squared statistics
		chiSquared <- apply(imputationGridResults, 2, function(x)
		{
			x[x > founders] <- NA
			subsettedX <- x[!finalSubsetIsNA]
			observed <- table(subsettedX, finalSubsetOmitNA)
			marginalImputed <- .rowSums(observed, m = nrow(observed), n = ncol(observed))
			marginalNewMarker <- .colSums(observed, m = nrow(observed), n = ncol(observed))
			expected <- outer(marginalImputed, marginalNewMarker) / sum(marginalImputed)
			return(sum((observed - expected)^2 / expected, na.rm=TRUE))
		})
		if(onlyStatistics)
		{
			testMap <- imputationMap
			testMap <- lapply(testMap, function(x) x[names(x) %in% extraImputationPoints(mpcrossMapped)])
			class(testMap) <- "map"
			return(new("addExtraMarkersStatistics", data = chiSquared, map = testMap))
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

		#Combine and keep recombination fractions
		combined <- combineKeepRF(mpcrossMapped, newMarkers, verbose = verbose, gbLimit = mpcrossMapped@rf@gbLimit, callEstimateRF = match(names(mpcrossMapped@map[[bestChromosome]]), markers(mpcrossMapped)), skipValidity = TRUE)
		#Put the new subset of markers in the right place.
		splitResults <- splitVector(names(relevantChromosomeMap)[markerRange], names(relevantChromosomeMap)[markerIndex])
		relevantSubset <- subset(combined, markers = c(splitResults$before, markers(newMarkers), splitResults$after), skipValidity = TRUE)
		grouped <- formGroups(relevantSubset, groups = 1, clusterBy = "theta")
		if(reorder)
		{
			if(verbose) cat("Reordering markers [", min(markerRange), ":", max(markerRange), "] of the chromosome, total number of markers for this chromosome was ", length(relevantChromosomeMap), "\n", sep = "")
#			if(hasMpMapInteractive2)
#			{
#				reorderedGrouped <- mpMapInteractive2::mpMapInteractive2(grouped)$object
#			}
#			else
#			{
				reorderedGrouped <- do.call(orderCross, c(list(grouped), orderCrossArgs))
#			}
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
		objectInNewOrder <- subset(combined, markers = newMarkerOrder, skipValidity = TRUE)
		newChromosomeMarkers <- newMarkerOrder[(newMarkerOrder %in% names(mpcrossMapped@map[[bestChromosome]])) | (newMarkerOrder %in% markers(newMarkers))]
		#Get out just the chromosome that's been changed. We only need to recompute the map for these markers.
		changedChromosomeInNewOrder <- subset(objectInNewOrder, markers = newChromosomeMarkers, skipValidity = TRUE)

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
		objectInNewOrder <- new("mpcrossMapped", objectInNewOrder, map = finalMap, rf = objectInNewOrder@rf, skipValidity = TRUE)

		#Update the imputation data, so we don't have to run the entire thing all over again. 
		if(!is.null(imputationArgs))
		{
			previousKey <- mpcrossMapped@geneticData[[1]]@imputed@key
			mappedChangedChromosomeInNewOrder <- new("mpcrossMapped", changedChromosomeInNewOrder, rf = changedChromosomeInNewOrder@rf, map = finalMap[bestChromosome], skipValidity = TRUE)
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
