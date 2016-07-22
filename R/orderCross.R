#' @export
orderCross <- function(mpcrossLG, cool = 0.5, tmin = 0.1, nReps = 1, maxMove = 0, effortMultiplier = 1, randomStart = TRUE, verbose = FALSE)
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
	permutation <- .Call("order", mpcrossLG, mpcrossLG@lg@allGroups, cool, tmin, nReps, maxMove, effortMultiplier, randomStart, verbose, PACKAGE="mpMap2")
	return(subset(mpcrossLG, markers = permutation))
}
#' @export
clusterOrderCross <- function(mpcrossLG, cool = 0.5, tmin = 0.1, nReps = 1, maxMove = 0, effortMultiplier = 1, randomStart = TRUE, nGroups)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	if(is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input mpcrossLG object did not contain recombination fraction information")
	}
	if(missing(nGroups) || nGroups < 1)
	{
		stop("Input nGroups must be a positive integer")
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
		if(nGroups >= length(underlying@markers))
		{
			orderedMarkers[[groupAsCharacter]] <- underlying@markers
		}
		else
		{
			distMatrix <- .Call("rawSymmetricMatrixToDist", underlying, PACKAGE="mpMap2")
			clustered <- fastcluster::hclust(distMatrix, method = "average")
			groupings <- cutree(clustered, nGroups)
			dissimilarityMatrix <- .Call("constructDissimilarityMatrix", underlying, groupings, PACKAGE="mpMap2")
	
			diag(dissimilarityMatrix) <- 0
			orderingOfGroups <- .Call("arsa", nGroups, dissimilarityMatrix[upper.tri(dissimilarityMatrix, diag = TRUE)], cool, tmin, nReps, maxMove, effortMultiplier, randomStart, PACKAGE="mpMap2")
			ordering <- unlist(lapply(1:orderingOfGroups, function(x) which(groupings == orderingOfGroups[x]+1)))
			orderedMarkers[[groupAsCharacter]] <- underlying@markers[ordering]
		}
	}
	return(subset(mpcrossLG, markers = unlist(orderedMarkers)))
}
