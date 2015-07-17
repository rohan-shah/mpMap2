setMethod(f = "subset", signature = "mpcross", definition = function(x, ...)
{
	newGeneticData <- lapply(x@geneticData, 
		function(geneticData, ...)
		{
			subset(geneticData, ...)
		}
	, ...)
	return(new("mpcross", geneticData = newGeneticData))
})
setMethod(f = "subset", signature = "mpcrossMapped", definition = function(x, ...)
{
	subsettedRF <- NULL
	if(!is.null(x@rf)) subsettedRF <- subset(x@rf, ...)
	return(new("mpcrossMapped", callNextMethod(), map = subset(x@map, ...), rf = subsettedRF))
})
setMethod(f = "subset", signature = "mpcrossRF", definition = function(x, ...)
{
	return(new("mpcrossRF", callNextMethod(), "rf" = subset(x@rf, ...)))
})
setMethod(f = "subset", signature = "mpcrossLG", definition = function(x, ...)
{
	subsettedRF <- NULL
	if(!is.null(x@rf)) subsettedRF <- subset(x@rf, ...)
	return(new("mpcrossLG", callNextMethod(), "lg" = subset(x@lg, ...), "rf" = subsettedRF))
})
setMethod(f = "subset", signature = "lg", definition = function(x, ...)
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
		markerIndices <- match(markers, markers(x))
	}

	groups <- x@groups[markerIndices]
	allGroups <- sort(unique(groups))
	return(new("lg", groups = groups, allGroups = allGroups))
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
		markerIndices <- match(markers, markers(x))
	}

	return(new("geneticData", founders = x@founders[,markerIndices,drop=FALSE], finals = x@finals[,markerIndices,drop=FALSE], hetData = subset(x@hetData, ...), pedigree = x@pedigree))
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
		markerIndices <- match(markers, markers(x))
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
		markerIndices <- match(markers, markers(x))
	}
	
	newTheta <- x@theta[markerIndices, markerIndices,drop=FALSE]
	
	if(is.null(x@lod))
	{
		newLod <- NULL
	}
	else
	{
		newLod <- x@lod[markerIndices, markerIndices,drop=FALSE]
	}

	if(is.null(x@lkhd))
	{
		newLkhd <- NULL
	}
	else
	{
		newLkhd <- x@lkhd[markerIndices, markerIndices,drop=FALSE]
	}
	return(new("rf", r = x@r, theta = newTheta, lod = newLod, lkhd = newLkhd))
})