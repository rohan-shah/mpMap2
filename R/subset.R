setMethod(f = "subset", signature = "mpcross", definition = function(x, ...)
{
		newGeneticData <- lapply(x@geneticData, 
			function(geneticData, ...)
			{
				subset(geneticData, ...)
			}
		, ...)
	newMap <- NULL
	if(!is.null(x@map))
	{
		newMap <- subset(x@map, ...)
	}
	return(new("mpcross", geneticData = newGeneticData, map = newMap))
})
setMethod(f = "subset", signature = "mpcrossRF", definition = function(x, ...)
{
	return(new("mpcrossRF", callNextMethod(), "rf" = subset(x@rf, ...)))
})
setMethod(f = "subset", signature = "geneticData", definition = function(x, ...)
{
	arguments <- list(...)
	if(!("markers" %in% names(arguments)))
	{
		stop("Argument markers is required for function subset.geneticData")
	}
	markers <- arguments$markers
	if(mode(markers) == "numeric")
	{
		markerIndices <- markers
	}
	else if(mode(markers) == "character")
	{
		markerIndices <- match(markers, colnames(x@founders))
	}

	return(new("geneticData", founders = x@founders[,markerIndices], finals = x@finals[,markerIndices], hetData = subset(x@hetData, ...), pedigree = x@pedigree))
})
setMethod(f = "subset", signature = "hetData", definition = function(x, ...)
{
	arguments <- list(...)
	if(!("markers" %in% names(arguments)))
	{
		stop("Argument markers is required for function subset.geneticData")
	}
	markers <- arguments$markers
	if(mode(markers) == "numeric")
	{
		markerIndices <- markers
	}
	else if(mode(markers) == "character")
	{
		markerIndices <- match(markers, colnames(x@founders))
	}
	return(new("hetData", x[markerIndices]))
})
setMethod(f = "subset", signature = "rf", definition = function(x, ...)
{
	arguments <- list(...)
	if(!("markers" %in% names(arguments)))
	{
		stop("Argument markers is required for function subset.rf")
	}
	markers <- arguments$markers
	if(mode(markers) == "numeric")
	{
		markerIndices <- markers
	}
	else if(mode(markers) == "character")
	{
		markerIndices <- match(markers, colnames(x@theta))
	}
	
	newTheta <- x@theta[markerIndices, markerIndices]
	
	if(is.null(x@lod))
	{
		newLod <- NULL
	}
	else
	{
		newLod <- x@lod[markerIndices, markerIndices]
	}

	if(is.null(x@lkhd))
	{
		newLkhd <- NULL
	}
	else
	{
		newLkhd <- x@lkhd[markerIndices, markerIndices]
	}
	return(new("rf", r = x@r, theta = newTheta, lod = newLod, lkhd = newLkhd))
})