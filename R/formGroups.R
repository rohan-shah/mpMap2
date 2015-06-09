formGroups <- function(mpcrossRF, groups, clusterBy="combined", method="average")
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
	nonNegativeIntegerArgument(groups)

	if((clusterBy %in% c("combined", "lod")) && is.null(mpcrossRF@rf@lod))
	{
		stop("Input mpcrossRF object must have a @rf@lod entry (likelihood ratio) in order to use combined or lod grouping")
	}
	if(clusterBy %in% c("combined", "lod"))
	{
		lod <- mpcrossRF@rf@lod
		#Reverse lod so that small values indicate similarity
		lod[is.na(mpcrossRF@rf@lod)] <- 0
		lod <- max(lod) - lod
		diag(lod) <- 0
	}
	theta <- mpcrossRF@rf@theta
	
	theta[is.na(mpcrossRF@rf@theta)] <- 0.5

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
		distMat <- theta + lod / max(lod) * min(abs(diff(mpcrossRF@rf@r)))
	}
	else if(clusterBy == "theta")
	{
		distMat <- theta
	}
	else
	{
		distMat <- lod
	}
	clustered <- hclust(as.dist(distMat), method=method)

	cut <- cutree(clustered, k=groups)
	names(cut) <- markers(mpcrossRF)
	
	lg <- new("lg", allGroups=1:groups, groups=cut)
	output <- new("mpcrossLG", mpcrossRF, lg = lg)
	return(subset(output, markers = order(cut)))
}