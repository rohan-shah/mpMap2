#' Estimate recombination fractions
#' 
#' This function estimates the recombination fractions between all pairs of markers in the input object. The recombination fractions are estimated using numerical maximum likelihood, and a grid search. Because every estimate will be one of the input test values, the estimates can be stored efficiently with a single byte per estimate.
#' @param object The input mpcross object
#' @param recombValues a vector of test values to use for the numeric maximum likelihood step. Must contain 0 and 0.5, and must have less than 255 values in total. The default value is \code{c(0:20/200, 11:50/100)}. 
#' @param lineWeights Values to use to correct for segregation distortion. This parameter should in general be left unspecified. 
#' @param gbLimit The maximum amount of working memory this estimation step should be allowed to use at any one time, in gigabytes. Smaller values may increase the computation time. A value of -1 indicates no limit.  
#' @param keepLod Set to \code{TRUE} to compute the likelihood ratio score statistics for testing whether the estimate is different from 0.5. Due to memory constraints this should generally be left as \code{FALSE}. 
#' @param keepLkhd Set to \code{TRUE} to compute the maximum value of the likelihood. Due to memory constraints this should generally be left as \code{FALSE}.
#' @param verbose Output diagnostic information, such as the amount of memory required, and the progress of the computation
#' @export
#' @examples map <- qtl::sim.map(len = 100, n.mar = 11, include.x=FALSE)
#' f2Pedigree <- f2Pedigree(1000)
#' cross <- simulateMPCross(map = map, pedigree = f2Pedigree, mapFunction = haldane, seed = 1)
#' rf <- estimateRF(cross)
#' #Print the estimated recombination fraction values
#' rf@@rf@@theta[1:11, 1:11]
estimateRF <- function(object, recombValues, lineWeights, gbLimit = -1, keepLod = FALSE, keepLkhd = FALSE, verbose = FALSE)
{
	inheritsNewMpcrossArgument(object)

	if (missing(recombValues)) recombValues <- c(0:20/200, 11:50/100)
	if (length(recombValues) >= 255)
	{
		stop("This package currently allows a maximum of 254 possible recombination fraction values")
	}
	recombValues <- sort(recombValues)
	if(!(0.5 %in% recombValues))
	{
		stop("Input recombValues must contain a value of 0.5")
	}
	if(!(0 %in% recombValues))
	{
		stop("Input recombValues must contain a value of 0")
	}
	if(missing(lineWeights))
	{
		lineWeights <- lapply(object@geneticData, function(x) rep(1, nLines(x)))
	}
	if(is.logical(verbose))
	{
		if(is.na(verbose))
		{
			stop("Input verbose cannot be NA")
		}
		else if(verbose)
		{
			verbose <- list(verbose = TRUE, progressStyle = 3L)
		}
		else verbose <- list(verbose = FALSE, progressStyle = 3L)
	}
	else
	{
		if(!is.list(verbose) || !("progressStyle" %in% names(verbose)) || !("verbose" %in% names(verbose)))
		{
			stop("Input verbose must be TRUE, FALSE, or a list with entries named progressStyle and verbose")
		}
		if(length(verbose$progressStyle) != 1L || !(verbose$progressStyle %in% 0:3))
		{
			stop("Input verbose$progressStyle must have value 0, 1, 2 or 3")
		}
		if(!is.logical(verbose$verbose) || length(verbose$verbose) != 1L)
		{
			stop("Input verbose$verbose must have value 0, 1, 2 or 3")
		}
	}

	if(class(lineWeights) == "numeric") lineWeights <- list(lineWeights)
	isNumericVectorListArgument(lineWeights)
	for(i in 1:length(object@geneticData))
	{
		if(length(lineWeights[[i]]) != nLines(object@geneticData[[i]]))
		{
			stop(paste0("Value of lineWeights[[", i, "]] must have nLines(object)[", i, "] entries"))
		}
	}
	markerRange <- 1:nMarkers(object)
	listOfResults <- estimateRFInternal(object = object, recombValues = recombValues, lineWeights = lineWeights, markerRows = markerRange, markerColumns = markerRange, keepLod = keepLod, keepLkhd = keepLkhd, gbLimit = gbLimit, verbose = verbose)
	theta <- new("rawSymmetricMatrix", markers = markers(object), levels = recombValues, data = listOfResults$theta)
	if(!is.null(listOfResults$lod))
	{
		listOfResults$lod <- new("dspMatrix", Dim = c(length(markers(object)), length(markers(object))), x = listOfResults$lod)
		rownames(listOfResults$lod) <- colnames(listOfResults$lod) <- markers(object)
	}
	if(!is.null(listOfResults$lkhd))
	{
		listOfResults$lkhd <- new("dspMatrix", Dim = c(length(markers(object)), length(markers(object))), x = listOfResults$lkhd)
		rownames(listOfResults$lkhd) <- colnames(listOfResults$lkhd) <- markers(object)
	}
	rf <- new("rf", theta = theta, lod = listOfResults$lod, lkhd = listOfResults$lkhd, gbLimit = gbLimit)

	if(class(object) == "mpcrossLG" || class(object) == "mpcrossMapped")
	{
		output <- object
		output@rf <- rf
	}
	else
	{
		output <- new("mpcrossRF", geneticData = object@geneticData, rf = rf)
	}
	return(output)
}
estimateRFInternal <- function(object, recombValues, lineWeights, markerRows, markerColumns, keepLod, keepLkhd, gbLimit, verbose)
{
	return(.Call("estimateRF", object, recombValues, markerRows, markerColumns, lineWeights, keepLod, keepLkhd, gbLimit, verbose, PACKAGE="mpMap2"))
}
