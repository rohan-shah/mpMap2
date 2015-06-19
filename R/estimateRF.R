estimateRF <- function(object, recombValues, lineWeights, keepLod = FALSE, keepLkhd = FALSE)
{
	isNewMpcrossArgument(object)

	if (missing(recombValues)) recombValues <- c(0:20/200, 11:50/100)
	if(!(0.5 %in% recombValues))
	{
		stop("Input recombValues must contain a value of 0.5")
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
	return(estimateRFSubset(object = object, recombValues=recombValues, lineWeights = lineWeights, start1 = 1, finish1 = nMarkers(object), start2 = 1, finish2 = nMarkers(object), keepLod = keepLod, keepLkhd = keepLkhd))
}
estimateRFInternal <- function(object, recombValues, lineWeights, marker1Range, marker2Range, keepLod, keepLkhd)
{
	return(.Call("estimateRF", object, recombValues, marker1Range, marker2Range, lineWeights, keepLod, keepLkhd, PACKAGE="mpMap2"))
}
estimateRFSubset <- function(object, recombValues, lineWeights, start1, finish1, start2, finish2, keepLod, keepLkhd)
{ 
	marker1Range <- c(start1,finish1)
	marker2Range <- c(start2,finish2)
	rpairs <- estimateRFInternal(object,  recombValues, lineWeights, marker1Range, marker2Range, keepLod, keepLkhd)
	rf <- new("rf", r = rpairs$r, theta = rpairs$theta, lod = rpairs$lod, lkhd = rpairs$lkhd)

	output <- new ("mpcrossRF", geneticData = object@geneticData, rf = rf)
	return(output)
}

