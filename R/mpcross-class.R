#' @include rf-class.R
#' @include lg-class.R
#' @include map-class.R
#' @include geneticData-class.R
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
	allLineNames <- unlist(lapply(object@geneticData, function(x) rownames(x@finals)))
	if(length(allLineNames) != length(unique(allLineNames)))
	{
		errors <- c(errors, "Line names must be unique")
	}
	#Check that the line names are unique
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.mpcross <- setClass("mpcross", slots = list(geneticData = "geneticDataList"), validity=checkMpcross)

checkMpcrossRF <- function(object)
{
	errors <- c()
	if(length(object@rf@theta@markers) != length(markers(object)))
	{
		errors <- "Inconsistent number of markers in slots @geneticData[[1]] and @rf@theta"
	}
	else if(any(object@rf@theta@markers != markers(object)))
	{
		errors <- "Inconsistent markers in slots @geneticData[[1]] and @rf@theta"
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
setClassUnion("rfOrNULL", c("rf", "NULL"))
.mpcrossRF <- setClass("mpcrossRF", contains = "mpcross", slots = list(rf = "rf"), validity=checkMpcrossRF)

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
	if(!is.null(object@lg@imputedTheta) && !is.null(object@rf))
	{
		if(length(object@lg@imputedTheta) > 0)
		{
			if(!identical(object@lg@imputedTheta[[1]]@levels, object@rf@theta@levels))
			{
				return("Slots lg@imputedTheta and rf@theta must have the same levels")
			}
		}
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
				errors <- c(errors, "Marker names in the genetic data must be the same as those in the map, and must occur in the same order as unlist(lapply(map, names))")
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
setAs("mpcrossLG", "mpcrossRF", def = function(from, to)
	{
		if(is.null(from@rf) && !is.null(from@lg@imputedTheta) && length(from@lg@allGroups) == 1)
		{
			return(new(to, as(from, "mpcross"), rf = new("rf", theta = from@lg@imputedTheta[[1]])))
		}
		else if(is.null(from@rf))
		{
			stop("As no RF data is present, this object of class mpcrossLG cannot be automatically converted to an object of class mpcrossRF. Please call estimateRF to re-estimate recombination fractions")
		}
		else return(new(to, as(from, "mpcross"), rf = from@rf))
	})
#' @export
infiniteSelfing <- function(founders, finals, pedigree)
{
	hetData <- lapply(1:ncol(founders), function(x)
		{
			alleles <- unique(founders[,x])
			cbind(alleles,alleles,alleles)
		})
	names(hetData) <- colnames(founders)
	hetData <- new("hetData", hetData)
	return(hetData)
}
#' @export
mpcross <- function(founders, finals, pedigree, hetData = infiniteSelfing, fixCodingErrors = FALSE)
{
	if(!isS4(pedigree) || !inherits(pedigree, "pedigree"))
	{
		stop("Input pedigree must be an S4 object of class peigree")
	}
	if(is.character(founders) || is.null(dim(founders)) || length(dim(founders)) != 2)
	{
		stop("Input founders must be a numeric matrix")
	}
	if(is.data.frame(founders)) founders <- as.matrix(founders)
	mode(founders) <- "integer"
	if(is.character(finals) || is.null(dim(finals)) || length(dim(finals)) != 2)
	{
		stop("Input finals must be a numeric matrix")
	}
	if(is.data.frame(finals)) finals <- as.matrix(finals)
	if(is.null(rownames(finals)) || is.null(colnames(finals)) || is.null(rownames(founders)) || is.null(colnames(founders)))
	{
		stop("Inputs founders and finals must have row and column names")
	}
	mode(finals) <- "integer"
	if(is.function(hetData))
	{
		hetData <- hetData(founders, finals, pedigree)
	}
	else if(!isS4(hetData) || !inherits(hetData, "hetData"))
	{
		stop("Input hetData must be an object of class hetData")
	}
	if(ncol(founders) != ncol(finals) || ncol(finals) != length(hetData))
	{
		stop("Inputs hetData, founders and finals must have the same number of markers")
	}
	sortedFounderMarkers <- sort(colnames(founders))
	sortedFinalMarkers <- sort(colnames(finals))
	sortedHetDataMarkers <- sort(names(hetData))
	if(any(sortedFounderMarkers != sortedFinalMarkers) || any(sortedFinalMarkers != sortedHetDataMarkers))
	{
		stop("Inputs founders, finals and hetData must have the same markers")
	}
	#Standardise marker order, if required
	if(any(colnames(founders) != colnames(finals)) || any(colnames(finals) != names(hetData)))
	{
		finals <- finals[,sortedFounderMarkers]
		hetData <- hetData[sortedFounderMarkers]
	}
	codingErrors <- listCodingErrors(founders = founders, finals = finals, hetData = hetData)
	if(length(codingErrors$null))
	{
		hetData[codingErrors$null] <- list(matrix(0L, 0, 3))
		finals[,codingErrors$null] <- NA
		warning(paste0("Removing data for ", length(codingErrors$null), " markers, because these markers have NA founder alleles"))
	}

	if(fixCodingErrors)
	{
		uniqueMarkers <- unique(codingErrors$finals[,"Column"])
		finals[, uniqueMarkers] <- NA
		warning(paste0("Removing data for ", length(uniqueMarkers), " markers, because fixCodingErrors = TRUE was specified. For less aggressive removal, use listCodingErrors"))
	}
	geneticData <- new("geneticData", founders = founders, hetData = hetData, finals = finals, pedigree = pedigree)
	mpcross <- new("mpcross", geneticData = new("geneticDataList", list(geneticData)))
	return(mpcross)
}
#' @export
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
