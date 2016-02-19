#' @export
orderCross <- function(mpcrossLG, cool = 0.5, tmin = 0.1, nReps = 1, verbose = FALSE)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	if(is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input mpcrossLG object did not contain recombination fraction information")
	}
	if(nMarkers(mpcrossLG) == 1)
	{
		return(mpcrossLG)
	}
	mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	permutation <- .Call("order", mpcrossLG, mpcrossLG@lg@allGroups, cool, tmin, nReps, verbose, PACKAGE="mpMap2")
	return(subset(mpcrossLG, markers = permutation))
}
#' @export
orderLargeCross <- function(mpcrossLG, maxSize = 1000, cool = 0.5, tmin = 0.1, nReps = 1, verbose = FALSE)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	if(is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input mpcrossLG object did not contain recombination fraction information")
	}
	mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	orderedMarkers <- vector(mode = "list")
	for(group in mpcrossLG@lg@allGroups)
	{
		groupAsCharacter <- as.character(group)
		currentGroupCross <- subset(mpcrossLG, groups = group)
		currentGroupCross <- subset(currentGroupCross, markers = sample(nMarkers(currentGroupCross)))
		if(verbose) cat("Started ordering group ", group, " which has ", nMarkers(currentGroupCross), "markers\n", sep="")
		if(nMarkers(currentGroupCross) < maxSize)
		{
			ordered <- orderCross(currentGroupCross)
			orderedMarkers[[groupAsCharacter]] <- markers(ordered)
			if(verbose) cat("This group had few markers\n", sep="")
			next
		}

		nGroups <- round(nMarkers(currentGroupCross) / maxSize)
		nGroups <- max(nGroups, 1)
		groupSize <- ceiling(nMarkers(currentGroupCross) / nGroups)
		starts <- round(seq(1, nMarkers(currentGroupCross), groupSize))
		if(nMarkers(currentGroupCross) %in% starts)
		{
			starts <- sort(setdiff(starts, nMarkers(currentGroupCross)))
		}
		ends <- c(starts[-1]-1, nMarkers(currentGroupCross))
		interleavedMarkers <- c()
		for(i in 1:length(starts))
		{
			part <- subset(currentGroupCross, markers = starts[i]:ends[i])
			if(verbose) cat("Initial ordering step ", i, " / ", length(starts), "\n", sep="")
			orderedPart <- orderCross(part, verbose = FALSE)
			if(length(interleavedMarkers) == 0)
			{
				interleavedMarkers <- markers(orderedPart)
			}
			else
			{
				newInterleaved <- vector(mode = "character", length = length(interleavedMarkers) + nMarkers(orderedPart))
				insertionIndices <- round(seq(1, length(newInterleaved), length.out = nMarkers(orderedPart)))
				newInterleaved[-insertionIndices] <- interleavedMarkers
				newInterleaved[insertionIndices] <- markers(orderedPart)

				adjacentIndices <- insertionIndices+1
				copiedInsertionIndices <- insertionIndices
				#There should only be one adjacent index (maybe) that has to be adjusted. 
				if(sum(adjacentIndices > length(newInterleaved)) > 0)
				{
					adjacentIndices <- adjacentIndices[-length(adjacentIndices)]
					copiedInsertionIndices <- copiedInsertionIndices[-length(copiedInsertionIndices)]
				}
				#consider reversing the order of the markers in the ordered part. 
				if(is.null(mpcrossLG@lg@imputedTheta))
				{
					adjacentIndices <- match(newInterleaved[adjacentIndices], mpcrossLG@rf@theta@markers)
					copiedInsertionIndices <- match(newInterleaved[copiedInsertionIndices], mpcrossLG@rf@theta@markers)
					currentOrderSum <- sum(mpcrossLG@rf@theta[cbind(adjacentIndices, copiedInsertionIndices)])
					oppositeOrderSum <- sum(mpcrossLG@rf@theta[cbind(adjacentIndices, rev(copiedInsertionIndices))])
				}
				else 
				{
					adjacentIndices <- match(newInterleaved[adjacentIndices], mpcrossLG@lg@imputedTheta[[as.character(group)]]@markers)
					copiedInsertionIndices <- match(newInterleaved[copiedInsertionIndices], mpcrossLG@lg@imputedTheta[[as.character(group)]]@markers)
					currentOrderSum <- sum(mpcrossLG@lg@imputedTheta[[as.character(group)]][cbind(adjacentIndices, copiedInsertionIndices)])
					oppositeOrderSum <- sum(mpcrossLG@lg@imputedTheta[[as.character(group)]][cbind(adjacentIndices, rev(copiedInsertionIndices))])
				}
				if(oppositeOrderSum < currentOrderSum) newInterleaved[insertionIndices] <- rev(newInterleaved[insertionIndices])
				interleavedMarkers <- newInterleaved
			}
		}
		orderedMarkers[[groupAsCharacter]] <- interleavedMarkers
		#So we've done the interleaving, now we need to order individual parts
		starts <- round(seq(1, nMarkers(currentGroupCross), maxSize / 2))
		starts <- starts[starts < nMarkers(currentGroupCross) - maxSize/2]
		for(start in starts)
		{
			end <- min(start + maxSize, nMarkers(currentGroupCross))
			originalMarkersThisChunk <- orderedMarkers[[groupAsCharacter]][start:end]
			object <- subset(currentGroupCross, markers = originalMarkersThisChunk)
			if(verbose) cat("Fine ordering step ", start, " / ", length(starts), "\n", sep="")
			ordered <- orderCross(object, verbose=FALSE)
			markersThisChunk <- markers(ordered)
			if(cor(match(markersThisChunk, originalMarkersThisChunk), 1:length(originalMarkersThisChunk)) < 0) markersThisChunk <- rev(markersThisChunk)
			#Again, make sure we don't accidentally reverse orderings. 
			orderedMarkers[[groupAsCharacter]][start:end] <- markersThisChunk
		}
	}
	return(subset(mpcrossLG, markers = unlist(orderedMarkers)))
}
