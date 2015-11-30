#' @export
estimateRF <- function(object, recombValues, lineWeights, keepLod = FALSE, keepLkhd = FALSE)
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
	listOfResults <- estimateRFInternal(object = object, recombValues = recombValues, lineWeights = lineWeights, markerRows = markerRange, markerColumns = markerRange, keepLod = keepLod, keepLkhd = keepLkhd)
	theta <- new("rawSymmetricMatrix", markers = markers(object), levels = recombValues, data = listOfResults$theta)
	if(!is.null(listOfResults$lod))
	{
		listOfResults$lod <- new("dspMatrix", Dim = c(length(markers(object)), length(markers(object))), x = listOfResults$lod)
		colnames(listOfResults$lod) <- markers(object)
	}
	if(!is.null(listOfResults$lkhd))
	{
		listOfResults$lkhd <- new("dspMatrix", Dim = c(length(markers(object)), length(markers(object))), x = listOfResults$lkhd)
		colnames(listOfResults$lkhd) <- markers(object)
	}
	rf <- new("rf", theta = theta, lod = listOfResults$lod, lkhd = listOfResults$lkhd)

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
estimateRFInternal <- function(object, recombValues, lineWeights, markerRows, markerColumns, keepLod, keepLkhd)
{
	return(.Call("estimateRF", object, recombValues, markerRows, markerColumns, lineWeights, keepLod, keepLkhd, PACKAGE="mpMap2"))
}
