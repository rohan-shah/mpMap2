#' @export
formGroups <- function(mpcrossRF, groups, clusterBy="combined", method="average", preCluster = FALSE)
{
	if(!(clusterBy %in% c("combined", "theta", "lod")))
	{
		stop("Input clusterBy must be one of 'combined', 'theta' or 'lod'")
	}
	if(!(method %in% c("average", "complete", "single")))
	{
		stop("Input method must be one of 'average', 'complete' or 'single'")
	}
	isNewMpcrossRFArgument(mpcrossRF)
	mpcrossRF <- as(mpcrossRF, "mpcrossRF")
	nonNegativeIntegerArgument(groups)

	if((clusterBy %in% c("combined", "lod")) && is.null(mpcrossRF@rf@lod))
	{
		stop("Input mpcrossRF object must have a @rf@lod entry (likelihood ratio) in order to use combined or lod grouping")
	}
	if(!preCluster)
	{
		nMarkers <- length(mpcrossRF@rf@theta@markers)
		theta <- mpcrossRF@rf@theta[1:nMarkers, 1:nMarkers]
		theta[is.na(theta)] <- 0.5
		if(clusterBy %in% c("combined", "lod"))
		{
			if(is.null(mpcrossRF@rf@lod))
			{
				stop("Lod must be calculated if clusterBy is \"combined\" or \"lod\"")
			}
			lod <- as(mpcrossRF@rf@lod, "matrix")
			#Reverse lod so that small values indicate similarity
			lod[is.na(lod)] <- 0
			lod <- max(lod) - lod
			diag(lod) <- 0
		}
		if(method == "average")
		{
			linkFunction <- function(x) mean(x, na.rm=TRUE)
		}
		else if(method == "complete")
		{
			linkFunction <- function(x) max(x, na.rm=TRUE)
		}
		else
		{
			linkFunction <- function(x) min(x, na.rm=TRUE)
		}
		if(clusterBy == "combined")
		{
			#Cluster by theta first, and then by lod to break any ties
			distMat <- as(lod / max(lod) * min(abs(diff(mpcrossRF@rf@theta@levels))), "matrix") + theta
		}
		else if(clusterBy == "theta")
		{
			distMat <- theta
		}
		else
		{
			distMat <- lod
		}
		clustered <- fastcluster::hclust(as.dist(distMat), method=method)
		cut <- cutree(clustered, k=groups)
		names(cut) <- markers(mpcrossRF)
	}
	#If we have a huge number of markers, it might be necessary to do a pre-clustering step, where we join together all the markers that have zero recombination fractions. 
	else
	{
		preClusterResults <- .Call("preClusterStep", mpcrossRF, PACKAGE="mpMap2")
		if(clusterBy == "combined")
		{
			distMat <- .Call("hclustCombinedMatrix", mpcrossRF, preClusterResults, PACKAGE="mpMap2")
		}
		else if(clusterBy == "theta")
		{
			distMat <- .Call("hclustThetaMatrix", mpcrossRF, preClusterResults, PACKAGE="mpMap2")
		}
		else
		{
			distMat <- .Call("hclustLodMatrix", mpcrossRF, preClusterResults, PACKAGE="mpMap2")
		}
		attr(distMat, "Size") <- length(preClusterResults)
		clustered <- fastcluster::hclust(distMat, method = method)
		#This cut is for the grouped markers, we want a similar object for the ungrouped markers
		cut <- cutree(clustered, k=groups)

		originalCut <- vector(mode = "integer", length = nMarkers(mpcrossRF))
		names(originalCut) <- markers(mpcrossRF)
		for(i in 1:groups)
		{
			originalCut[unlist(preClusterResults[which(cut == i)])] <- i
		}
		cut <- originalCut
	}
	
	lg <- new("lg", allGroups=1:groups, groups=cut)
	output <- new("mpcrossLG", mpcrossRF, lg = lg, rf = mpcrossRF@rf)
	return(subset(output, markers = order(cut)))

}
