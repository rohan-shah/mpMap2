#' Computes map distances
#' 
#' Use the non-linear least squares function to estimate a map
#' @export
estimateMap <- function(mpcrossLG, mapFunction = rfToHaldane, maxOffset = 1)
{
	isNewMpcrossLGArgument(mpcrossLG)
	if (is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input object must have recombination fractions")
	}
	if(is.null(mpcrossLG@lg@imputedTheta))
	{
		cat("Imputing recombination fractions\n", sep="")
		mpcrossLG <- impute(mpcrossLG, verbose=TRUE)
	}
	map <- list()
	for (group in mpcrossLG@lg@allGroups)
	{
		rfData <- mpcrossLG@lg@imputedTheta[[as.character(group)]]
		maxOffset <- min(maxOffset,length(rfData@markers)-1)
		#Construct design matrix
		#d <- designMat(length(object$map[[chr]])-1, maxOffset)
		d <- .Call("generateDesignMatrix", length(rfData@markers)-1, maxOffset, PACKAGE="mpMap2")
		indices <- matrix(nrow=maxOffset * length(rfData@markers) - maxOffset * (maxOffset + 1) / 2, ncol = 2)
		counter <- 1
		for(offset in 1:maxOffset)
		{
			for(i in 1:(length(rfData@markers)-offset))
			{
				indices[counter,] <- c(i+offset, i)
				counter <- counter + 1
			}
		}
		#B vector for nnls
		b <- rfData[indices]
		#Values of 0.5 result in infinite estimated distance, which doesn't really work. 
		b[b == 0.5] <- 0.49
		if (!requireNamespace("nnls", quietly=TRUE))
			stop("nnls needed for computemap to work. Please install it.\n", call.=FALSE) 
		result <- nnls::nnls(d, mapFunction(b)) 
		map[[as.character(group)]] <- c(0, cumsum(result$x[which(indices[,1] == indices[,2]+1)]))
		names(map[[as.character(group)]]) <- rfData@markers
	}
	class(map) <- "map"
	return(map)
}
designMat <- function(n, maxOffset)
{
	resultMat <- matrix(0, nrow=n*maxOffset - maxOffset*(maxOffset - 1)/2, ncol=n)
	#column of matrix
	for(i in 1:n)
	{
		offset <- 1
		#j is the section going by rows
		#for(j in 1:n)
		for(j in 1:maxOffset)
		{
			resultMat[offset + max(0, i-j):min(n-j, i-1) ,i] <- 1
			offset <- offset + (n-j+1)
		}
	}
	return(resultMat)
}
