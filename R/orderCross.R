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
orderLargeCross <- function(mpcrossLG, cool = 0.5, tmin = 0.1, nReps = 1, markersPerGroup)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	if(is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input mpcrossLG object did not contain recombination fraction information")
	}
	if(missing(markersPerGroup) || markersPerGroup < 1)
	{
		stop("Input markersPerGroup must be a positive integer")
	}
	mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	orderedMarkers <- vector(mode = "list")
	for(group in mpcrossLG@lg@allGroups)
	{
		groupAsCharacter <- as.character(group)
		if(!is.null(mpcrossLG@lg@imputedTheta))
		{
			underlying <- mpcrossLG@lg@imputedTheta[[groupAsCharacter]]
		}
		else
		{
			currentGroupCross <- subset(mpcrossLG, groups = group)
			currentGroupCross <- impute(currentGroupCross)
			underlying <- currentGroupCross@lg@imputedTheta[[1]]
		}
		if(markersPerGroup >= length(underlying@markers))
		{
			orderedMarkers[[groupAsCharacter]] <- underlying@markers
		}
		else
		{
			distMatrix <- .Call("rawSymmetricMatrixToDist", underlying, PACKAGE="mpMap2")
			clustered <- fastcluster::hclust(distMatrix, method = "average")
			nGroups <- ceiling(length(underlying@markers) / markersPerGroup)
			groupings <- cutree(clustered, nGroups)
			dissimilarityMatrix <- .Call("constructDissimilarityMatrix", underlying, groupings, PACKAGE="mpMap2")
	
			diag(dissimilarityMatrix) <- 0
			orderingOfGroups <- .Call("arsa", nGroups, dissimilarityMatrix[upper.tri(dissimilarityMatrix, diag = TRUE)], cool, tmin, nReps, PACKAGE="mpMap2")
			ordering <- unlist(lapply(1:orderingOfGroups, function(x) which(groupings == orderingOfGroups[x]+1)))
			orderedMarkers[[groupAsCharacter]] <- underlying@markers[ordering]
		}
	}
	return(subset(mpcrossLG, markers = unlist(orderedMarkers)))
}
