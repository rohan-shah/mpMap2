#' @include rf-class.R
#' @include lg-class.R
#' @include map-class.R
checkCompatibleGeneticData <- function(objects)
{
	expectedMarkers <- markers(objects[[1]])
	errors <- c()
	for(i in 1:length(objects))
	{
		x <- objects[[i]]
		if(nMarkers(x) != length(expectedMarkers))
		{
			errors <- c(errors, paste0("Wrong number of markers in slot geneticData[[", i, "]]"))
		}
		else if(any(markers(x) != expectedMarkers))
		{
			errors <- c(errors, paste0("Wrong markers in slot geneticData[[", i, "]]"))
		}
	}

	return(errors)
}
checkMpcross <- function(object)
{
	errors <- c()
	if(length(object@geneticData) == 0)
	{
		errors <- c(errors, "Must contain at least one set of genetic data")
	}
	else
	{
		for(i in 1:length(object@geneticData))
		{
			x <- object@geneticData[[i]]
			if(class(x) != "geneticData")
			{
				errors <- c(errors, paste0("Value of slot geneticData[[", i, "]] must be a geneticData object"))
			}
			else
			{
				geneticErrors <- checkGeneticData(x)
				if(mode(geneticErrors) != "logical")
				{
					errors <- c(errors, paste0("Error in geneticData[[", i, "]]: ", geneticErrors))
				}
			}
		}
		errors <- c(errors, checkCompatibleGeneticData(object@geneticData))
	}
	#Check that the line names are unique
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.mpcross <- setClass("mpcross", slots = list(geneticData = "list"), validity=checkMpcross)

setClassUnion("rfOrNULL", c("rf", "NULL"))
.mpcrossRF <- setClass("mpcrossRF", contains = "mpcross", slots = list(rf = "rf"))

checkMpcrossLG <- function(object)
{
	if(is.null(names(object@lg@groups)))
	{
		return("Slot lg@groups must have names")
	}
	if(any(is.na(names(object@lg@groups)) || names(object@lg@groups) != markers(object)))
	{
		return("Marker names implied by names of slots lg@groups and founders were different")
	}
	return(TRUE)
}
.mpcrossLG <- setClass("mpcrossLG", contains = "mpcross", slots = list(lg = "lg", rf = "rfOrNULL"), validity=checkMpcrossLG)

checkMpcrossMapped <- function(object)
{
	errors <- c()
	nMapMarkers <- sum(unlist(lapply(object@map, length)))
	if(nMarkers(object) != nMapMarkers)
	{
		errors <- c(errors, "Number of markers in map is different from the number of markers in slot founders")
	}
	else
	{
		markerNames <- unlist(lapply(object@map, names))
		for(i in 1:length(object@geneticData))
		{
			geneticDataMarkers <- markers(object@geneticData[[i]])
			if(length(geneticDataMarkers) != length(markerNames) || any(geneticDataMarkers != markerNames))
			{
				errors <- c(errors, "Marker names in map disagree with marker names in genetic data")
			}
		}
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.mpcrossMapped <- setClass("mpcrossMapped", contains = "mpcross", slots = list(map = "map", rf = "rfOrNULL"), validity=checkMpcrossMapped)

setAs("mpcrossMapped", "mpcrossLG", def = function(from, to)
	{
		groups <- rep(1:length(from@map), each = unlist(lapply(from@map, length)))
		names(groups) <- markers(from)
		allGroups <- unique(groups)
		lg <- new("lg", allGroups = allGroups, groups = groups)
		return(new(to, as(from, "mpcross"), lg = lg, rf = from@rf))
	})
mpcrossMapped <- function(cross, map, rf=NULL)
{
	if(inherits(cross, "mpcrossRF"))
	{
		if(!is.null(rf))
		{
			stop("Two objects of class rf were specified")
		}
		return(new("mpcrossMapped", as(cross, "mpcross"), rf = cross@rf, map = map))
	}
	else
	{
		return(new("mpcrossMapped", cross, map = map, rf = rf))
	}
}
