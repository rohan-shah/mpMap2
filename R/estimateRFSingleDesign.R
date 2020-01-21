#' @export
#' @title Estimate pairwise recombination fractions
#' @description Estimate pairwise recombination fractions, similar to \code{\link{estimateRF}}, but with different performance requirements in terms of compute time and storage. 
#' @details Estimate pairwise recombination fractions, similar to \code{\link{estimateRF}}, but with different performance requirements in terms of compute time and storage. Specifically, this version is expected to perform better when there is only a single population. 
#' @param object An object of class \code{mpcross}.
#' @param recombValues a vector of test values to use for the numeric maximum likelihood step. Must contain 0 and 0.5, and must have less than 255 values in total. The default value is \code{c(0:20/200, 11:50/100)}.
#' @param lineWeights Values to use to correct for segregation distortion. This parameter should in general be left unspecified.
#' @param keepLod Set to \code{TRUE} to compute the likelihood ratio score statistics for testing whether the estimate is different from 0.5. Due to memory constraints this should generally be left as \code{FALSE}.
#' @param keepLkhd Set to \code{TRUE} to compute the maximum value of the likelihood. Due to memory constraints this should generally be left as \code{FALSE}.
#' @param verbose Output diagnostic information, such as the amount of memory required, and the progress of the computation.
#' @param markerRows Used to estimate only a subset of the full matrix of pairwise recombination fractions.
#' @param markerColumns Used to estimate only a subset of the full matrix of pairwise recombination fractions.
#' @return An object of class \code{mpcrossRF}, which contains the original genetic data, and also estimated recombination fraction data.
estimateRFSingleDesign <- function(object, recombValues, lineWeights, keepLod = FALSE, keepLkhd = FALSE, verbose = FALSE, markerRows = 1:nMarkers(object), markerColumns = 1:nMarkers(object))
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
		lineWeights <- rep(1, nLines(object))
	}
	else
	{
		stop("Use of lineWeights is not implented yet")
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
		if(!is.list(verbose) || !("progressStyle" %in% names(verbose)))
		{
			stop("Input verbose must be TRUE, FALSE, or a list with entries named progressStyle and verbose")
		}
		if(!("verbose" %in% names(verbose))) verbose$verbose <- TRUE
		if(length(verbose$progressStyle) != 1L || !(verbose$progressStyle %in% 0:3))
		{
			stop("Input verbose$progressStyle must have value 0, 1, 2 or 3")
		}
		if(!is.logical(verbose$verbose) || length(verbose$verbose) != 1L)
		{
			stop("Input verbose$verbose must have value 0, 1, 2 or 3")
		}
	}
	if(length(lineWeights) != nLines(object))
	{
		stop(paste0("Value of lineWeights  must have nLines(object) entries"))
	}
	listOfResults <- estimateRFSingleDesignInternal(object = object, recombValues = recombValues, lineWeights = lineWeights, markerRows = markerRows, markerColumns = markerColumns, keepLod = keepLod, keepLkhd = keepLkhd, verbose = verbose)
	if(is(object, "mpcrossRF") || ((is(object, "mpcrossLG") || is(object, "mpcrossMapped")) && !is.null(object@rf)))
	{
		warning("Updating RF data for existing object")
		.Call("assignRawSymmetricMatrixFromEstimateRF", object@rf@theta, markerRows, markerColumns, listOfResults$theta, PACKAGE="mpMap2")

		if(!is.null(listOfResults$lod))
		{
			if(!is.null(object@rf@lod)) lod <- object@rf@lod
			else lod <- new("dspMatrix", Dim = rep(nMarkers(object), 2), x = rep(as.numeric(NA), nMarkers(object)*(nMarkers(object) + 1)/2))
			.Call("assignDspMatrixFromEstimateRF", lod, markerRows, markerColumns, listOfResults$lod, PACKAGE="mpMap2")
			rownames(lod) <- colnames(lod) <- markers(object)
			object@rf@lod <- lod
		}

		if(!is.null(listOfResults$lkhd))
		{
			if(!is.null(object@rf@lkhd)) lkhd <- object@rf@lkhd
			else lkhd <- new("dspMatrix", Dim = rep(nMarkers(object), 2), x = rep(as.numeric(NA), nMarkers(object)*(nMarkers(object) + 1)/2))
			.Call("assignDspMatrixFromEstimateRF", lkhd, markerRows, markerColumns, listOfResults$lkhd, PACKAGE="mpMap2")
			rownames(lkhd) <- colnames(lkhd) <- markers(object)
			object@rf@lkhd <- lkhd
		}
		object@rf@gbLimit <- -1
		output <- object
	}
	else
	{
		theta <- new("rawSymmetricMatrix", markers = markers(object), levels = recombValues, data = rep(as.raw(0xFF), nMarkers(object)*(nMarkers(object) + 1)/2))
		.Call("assignRawSymmetricMatrixFromEstimateRF", theta, markerRows, markerColumns, listOfResults$theta, PACKAGE="mpMap2")

		lod <- NULL
		if(!is.null(listOfResults$lod))
		{
			lod <- new("dspMatrix", Dim = rep(nMarkers(object), 2), x = rep(as.numeric(NA), nMarkers(object)*(nMarkers(object) + 1)/2))
			.Call("assignDspMatrixFromEstimateRF", lod, markerRows, markerColumns, listOfResults$lod, PACKAGE="mpMap2")
			rownames(lod) <- colnames(lod) <- markers(object)
		}

		lkhd <- NULL
		if(!is.null(listOfResults$lkhd))
		{
			lkhd <- new("dspMatrix", Dim = rep(nMarkers(object), 2), x = rep(as.numeric(NA), nMarkers(object)*(nMarkers(object) + 1)/2))
			.Call("assignDspMatrixFromEstimateRF", lkhd, markerRows, markerColumns, listOfResults$lkhd, PACKAGE="mpMap2")
			rownames(lkhd) <- colnames(lkhd) <- markers(object)
		}
		rf <- new("rf", theta = theta, lod = lod, lkhd = lkhd, gbLimit = -1)

		if(class(object) == "mpcrossLG" || class(object) == "mpcrossMapped")
		{
			output <- object
			output@rf <- rf
		}
		else
		{
			output <- new("mpcrossRF", geneticData = object@geneticData, rf = rf)
		}
	}
	return(output)
}
estimateRFSingleDesignInternal <- function(object, recombValues, lineWeights, markerRows, markerColumns, keepLod, keepLkhd, verbose)
{
	return(.Call("estimateRFSingleDesign", object, recombValues, markerRows, markerColumns, lineWeights, keepLod, keepLkhd, verbose, PACKAGE="mpMap2"))
}
