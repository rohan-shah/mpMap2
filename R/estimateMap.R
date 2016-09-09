#' Computes map distances
#' 
#' @param maxMarkers The (approximate) number of markers for which distances are estimated simultaneously. 
#' Use the non-linear least squares function to estimate a map
#' @export
estimateMap <- function(mpcrossLG, mapFunction = rfToHaldane, maxOffset = 1, maxMarkers = 2000, verbose=FALSE)
{
	isNewMpcrossLGArgument(mpcrossLG)
	if (is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input object must have recombination fractions")
	}
	if(is.null(mpcrossLG@lg@imputedTheta))
	{
		cat("Imputing recombination fractions\n", sep="")
		mpcrossLG <- impute(mpcrossLG, verbose=verbose)
	}
	if(maxOffset > maxMarkers / 8)
	{
		stop("Input maxOffset must be smaller than maxMarkers / 8")
	}
	map <- list()
	for (group in mpcrossLG@lg@allGroups)
	{
		if(verbose) cat("Starting linkage group ", group, "\n", sep="")
		rfData <- mpcrossLG@lg@imputedTheta[[as.character(group)]]
		maxOffset <- min(maxOffset,length(rfData@markers)-1)
		originalMarkers <- length(rfData@markers)
		if(originalMarkers < maxMarkers)
		{
			ranges <- cbind(1, originalMarkers)
		}
		else
		{
			starts <- as.integer(seq(1, originalMarkers+1, length.out = ceiling(originalMarkers / maxMarkers)+1))
			ranges <- cbind(head(starts, -1), tail(starts, -1)-1)
		}
		resultsThisGroup <- list()
		if(verbose) cat("Splitting into ", length(ranges), " subproblems\n", sep="")
		for(row in 1:nrow(ranges))
		{
			start <- ranges[row, 1]
			end <- ranges[row, 2]
			if(verbose) 
			{
				cat("Iteration ", row, " / ", nrow(ranges), "\n", sep="")
			}
			subsettedRF <- subset(rfData, markers = start:end)
			resultsThisGroup[[row]] <- estimateMapInternal(subsettedRF = subsettedRF, maxOffset = maxOffset, verbose = verbose, mapFunction = mapFunction)
		}
		currentChr <- resultsThisGroup[[1]]
		#Now to glue all these bits together. 
		if(nrow(ranges) > 1)
		{
			for(row in 2:nrow(ranges))
			{
				#Get out the map for the last maxOffset+1 markers of the end of the last chunk, and the first maxOffset+1 markers of the next chunk. 
				subsettedRF <- subset(rfData, markers = (ranges[row-1,2] - 2*maxOffset):(ranges[row,1] + 2*maxOffset))
				overlap <- estimateMapInternal(subsettedRF, maxOffset, FALSE, mapFunction = mapFunction)
				currentChr[(length(currentChr) - maxOffset):(length(currentChr))] <- overlap[(maxOffset+1):(2*maxOffset+1)] - overlap[maxOffset+1] + currentChr[length(currentChr) - maxOffset]
				resultsThisGroup[[row]] <- resultsThisGroup[[row]] + tail(currentChr, 1) + (overlap[2*maxOffset+2] - overlap[2*maxOffset+1])
				currentChr <- c(currentChr, resultsThisGroup[[row]])
			}
		}
		map[[as.character(group)]] <- currentChr
	}
	class(map) <- "map"
	return(map)
}
estimateMapInternal <- function(subsettedRF, maxOffset, verbose, mapFunction)
{
	if(verbose) 
	{
		cat("Generating design matrx\n", sep="")
	}
	#Construct design matrix
	d <- .Call("generateDesignMatrix", length(subsettedRF@markers)-1, maxOffset, PACKAGE="mpMap2")
	indices <- matrix(nrow=maxOffset * length(subsettedRF@markers) - maxOffset * (maxOffset + 1) / 2, ncol = 2)
	counter <- 1
	for(offset in 1:maxOffset)
	{
		for(i in 1:(length(subsettedRF@markers)-offset))
		{
			indices[counter,] <- c(i+offset, i)
			counter <- counter + 1
		}
	}
	#B vector for nnls
	b <- subsettedRF[indices]
	#Values of 0.5 result in infinite estimated distance, which doesn't really work. 
	b[b == 0.5] <- 0.49
	if(verbose) cat("Solving non-negative least squares problem\n", sep="")
	result <- nnls::nnls(d, mapFunction(b))
	results <- c(0, cumsum(result$x[which(indices[,1] == indices[,2]+1)]))
	names(results) <- subsettedRF@markers
	return(results)
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
