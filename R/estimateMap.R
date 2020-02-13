#' @title Estimate map distances
#' @description Estimate map distances based on the estimated recombination fractions
#' @param mpcrossLG An object of class \code{mpcrossLG}, which must also contain data about recombination fractions and linkage groups. 
#' @param maxMarkers The (approximate) number of markers for which distances are estimated simultaneously. 
#' @param mapFunction The map function to use to compute recombination fractions to centiMorgan distances.
#' @param maxOffset The maximum separation between pairs of markers used for map construction, in terms of position within the ordering. Recombination fractions between pairs of markers, which are further apart than this, will not be used to estimate the map distances. 
#' @param verbose Should verbose output be produced?
#' @return A map object, in the format specified by the \code{\link[qtl]{qtl-package}} package. This format is a list of chromosomes, with each entry being a named numeric vector of marker positions. 
#' @details 
#' Once a marker order has been chosen, one possible way of estimating a genetic map is to convert the recombination fractions between adjacent markers into centiMorgan distances. This tends not to work well, because individual recombination fraction estimates can be highly variable, depending on the experimental design used, and the distribution of the marker alleles among the founders. It also wastes much of the information contained in the data; we can estimate recombination fractions between all pairs of markers, rather than just adjacent markers, and this information should be used in the estimation of map distances
#' 
#' This function uses non-linear least squares to estimate map distances as follows. Assume that there are \code{n} markers on a chromosome, and for all pairs of markers there is an available estimate of the recombination fraction. For every pair of markers which differ by \code{maxOffset} or less, in terms of their position within the ordering, the recombination fraction between these markers is turned into a centiMorgan distance. This centiMorgan distance is expressed as a sum of distances between adjacent markers, which is a simple equation. The set of all the equations generated in this way is represented as a matrix equation, and solved via non-linear least squares. As these non-linear least squares problems can become very large, input \code{maxMarkers} allows the non-linear least squares problem to be broken into several smaller problems. 
#'
#' For example, assume that there are five markers, for which an order has been determined. The distance between markers \eqn{i} and \eqn{j}, \emph{as estimated by the recombination fractions}, is \eqn{d(i, j)}. The genetic distance between markers \eqn{i} and \eqn{i + 1} \emph{in the final genetic map} is \eqn{a(i)}. So in this case, the parameters that are to be estimated are \eqn{a(1), a(2), a(3)} and \eqn{a(4)}. If \code{maxOffset} is 3, then the set of equations generated is 
#' \deqn{d(1, 3) = a(1) + a(2)}
#' \deqn{d(1, 4) = a(1) + a(2) + a(3)}
#' \deqn{d(2, 4) = a(2) + a(3)}
#' \deqn{d(3, 5) = a(3) + a(4)}
#' \deqn{d(2, 5) = a(2) + a(3) + a(4)}
#' These constraints are represented as a matrix equation and solved for \eqn{a(1), a(2), a(3)} and \eqn{a(4)} using non-linear least squares. However, if \code{maxOffset} is set to \code{2}, then the set of equations is 
#' \deqn{d(1, 3) = a(1) + a(2)}
#' \deqn{d(2, 4) = a(2) + a(3)}
#' \deqn{d(3, 5) = a(3) + a(4)}
#' @examples
#' data(simulatedFourParentData)
#' #Estimate recombination fractions
#' rf <- estimateRF(simulatedFourParentData)
#' #Assign all markers to one linkage group / chromosome
#' grouped <- formGroups(rf, groups = 1)
#' #Estimate map
#' \dontrun{estimatedMap <- estimateMap(grouped, maxOffset = 10)}
#' #Create object that includes the map
#' \dontrun{mapped <- new("mpcrossMapped", grouped, map = estimatedMap)}
#' @export
estimateMap <- function(mpcrossLG, mapFunction = rfToHaldane, maxOffset = 1, maxMarkers = 2000, verbose=FALSE)
{
	if(identical(mapFunction, haldane))
	{
		mapFunction <- rfToHaldane
	}
	if(identical(mapFunction, kosambi))
	{
		mapFunction <- rfToKosambi
	}
	isNewMpcrossLGArgument(mpcrossLG)
	if(!inherits(mpcrossLG, "mpcrossLG")) mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	if (is.null(mpcrossLG@rf) && is.null(mpcrossLG@lg@imputedTheta))
	{
		stop("Input object must have recombination fractions")
	}
	if(is.null(mpcrossLG@lg@imputedTheta))
	{
		if(verbose) cat("Imputing recombination fractions\n", sep="")
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
		if(verbose) cat("Splitting into ", nrow(ranges), " subproblems\n", sep="")
		for(row in 1:nrow(ranges))
		{
			start <- ranges[row, 1]
			end <- ranges[row, 2]
			if(verbose) 
			{
				cat("Iteration ", row, " / ", nrow(ranges), "\n", sep="")
			}
			subsettedRF <- subset(rfData, markers = start:end)
			resultsThisGroup[[row]] <- estimateMapInternal(subsettedRF = subsettedRF, maxOffset = maxOffset, mapFunction = mapFunction)
		}
		currentChr <- resultsThisGroup[[1]]
		#Now to glue all these bits together. 
		if(nrow(ranges) > 1)
		{
			for(row in 2:nrow(ranges))
			{
				#Get out the map for the last maxOffset+1 markers of the end of the last chunk, and the first maxOffset+1 markers of the next chunk. 
				subsettedRF <- subset(rfData, markers = (ranges[row-1,2] - 2*maxOffset):(ranges[row,1] + 2*maxOffset))
				overlap <- estimateMapInternal(subsettedRF = subsettedRF, maxOffset = maxOffset, mapFunction = mapFunction)
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
estimateMapInternal <- function(subsettedRF, maxOffset, mapFunction)
{
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
