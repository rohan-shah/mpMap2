checkMap <- function(object)
{
	allNumeric <- unlist(lapply(object, is.numeric))
	if(!allNumeric) return("A map object must be a list of numeric vectors")
	return(TRUE)
}
.map4 <- setClass("map4", contains = "list", validity = checkMap)
setOldClass("map", S4Class = "map4")
removeClass("map4")
setClassUnion("mapOrNULL", c("map", "NULL"))

checkPedigree <- function(object)
{
	nTotalLines <- length(object@lineNames)
	errors <- c()
	if(length(object@mother) != nTotalLines || length(object@father) != nTotalLines)
	{
		errors <- c(errors, "Lengths of slots lineNames, mother and father must be the same")
	}
	if(any(is.na(object@lineNames)))
	{
		errors <- c(errors, "Slot lineNames contained NA values")
	}
	if(any(is.na(object@mother)))
	{
		errors <- c(errors, "Slot mother contained NA values")
	}
	if(any(is.na(object@father)))
	{
		errors <- c(errors, "Slot father contained NA values")
	}

	if(any(object@mother < 0 | object@mother > nTotalLines))
	{
		errors <- c(errors, "Values in slot mother had invalid values")
	}
	if(any(object@father < 0 | object@father > nTotalLines))
	{
		errors <- c(errors, "Values in slot father had invalid values")
	}

	#Parents must come first in the pedigree
	if(!all(object@mother < 1:nTotalLines) || !all(object@father < 1:nTotalLines))
	{
		errors <- c(errors, "Mother and father must preceed offspring in the pedigree")
	}
	#Entry selfing must be either "auto" of "infinite"
	if(object@selfing != "auto" && object@selfing != "infinite")
	{
		errors <- c(errors, "Slot selfing must be either \"infinite\" or \"auto\"")
	}

	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.pedigree <- setClass("pedigree", slots = list(lineNames = "character", mother = "integer", father = "integer", selfing = "character"), validity = checkPedigree)
checkDetailedPedigree <- function(object)
{
	nTotalLines <- length(object@lineNames)
	errors <- c()
	if(length(object@observed) != nTotalLines)
	{
		errors <- c(errors, "Lengths of slots lineNames, mother, father and observed must be the same")
	}

	if(any(is.na(object@initial)))
	{
		errors <- c(errors, "Slot initial contained NA values")
	}
	if(any(is.na(object@observed)))
	{
		errors <- c(errors, "Slot observed contained NA values")
	}

	if(length(object@initial) == 0)
	{
		errors <- c(errors, "Slot initial must be non-empty")
	}

	if(any(object@initial < 0 | object@initial > nTotalLines))
	{
		errors <- c(errors, "Values in slot initial had invalid values")
	}

	if(any(object@mother[object@initial] != 0))
	{
		errors <- c(errors, "An entry in slot mother for a line named in slot initial cannot have a mother")
	}
	if(any(object@father[object@initial] != 0))
	{
		errors <- c(errors, "An entry in slot father for a line named in slot initial cannot have a father")
	}

	if(length(unique(object@initial)) != length(object@initial))
	{
		errors <- c(errors, "Slot initial cannot contain duplicate values")
	}

	if(!all(sort(object@initial) == 1:max(object@initial)))
	{
		errors <- c(errors, "Initial lines must be at the start of the pedigree")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.detailedPedigree <- setClass("detailedPedigree", contains = "pedigree", slots = list(initial = "integer", observed = "logical"), validity = checkDetailedPedigree)
checkHets <- function(object)
{
	if(is.null(names(object)) || any(names(object) == ""))
	{
		return("Every entry must have a valid name")
	}
	if(length(unique(names(object))) != length(object))
	{
		return("Names must be unique")
	}
	return(.Call("checkHets", object))
}
.hetData <- setClass("hetData", contains = "list", validity = checkHets)

checkGeneticData <- function(object)
{
	errors <- c()
	if(!is.integer(object@founders))
	{
		errors <- c(errors, "Slot founders must be an integer matrix")
	}
	if(!is.integer(object@finals))
	{
		errors <- c(errors, "Slot finals must be an integer matrix")
	}

	if(length(dim(object@finals)) != 2)
	{
		errors <- c(errors, "Slot finals must be a matrix")
	}
	if(length(dim(object@founders)) != 2)
	{
		errors <- c(errors, "Slot founders must be a matrix")
	}
	#Check that row and column names of finals and founders exist
	if(any(unlist(lapply(dimnames(object@founders), is.null))))
	{
		errors <- c(errors, "Slot founders must have row and column names")
	}
	if(any(unlist(lapply(dimnames(object@finals), is.null))))
	{
		errors <- c(errors, "Slot finals must have row and column names")
	}

	#Check for NA's in row and column names of finals and founders
	if(any(unlist(lapply(dimnames(object@founders), function(x) any(is.na(x))))))
	{
		errors <- c(errors, "Row and column names of slot founders cannot be NA")
	}
	if(any(unlist(lapply(dimnames(object@finals), function(x) any(is.na(x))))))
	{
		errors <- c(errors, "Row and column names of slot finals cannot be NA")
	}

	#Check that row and column names of finals and founders are unique
	if(any(unlist(lapply(dimnames(object@finals), function(x) length(x) != length(unique(x))))))
	{
		errors <- c(errors, "Row and column names of slot finals cannot contain duplicates")
	}
	if(any(unlist(lapply(dimnames(object@founders), function(x) length(x) != length(unique(x))))))
	{
		errors <- c(errors, "Row and column names of slot founders cannot contain duplicates")
	}

	nMarkers <- ncol(object@founders)
	markers <- colnames(object@finals)
	if(ncol(object@finals) != nMarkers || length(object@hetData) != nMarkers)
	{
		errors <- c(errors, "Slots finals, founders and hetData had different numbers of markers")
	}
	if(any(colnames(object@founders) != markers))
	{
		errors <- c(errors, "Slot finals must have the same colnames as slot founders")
	}
	if(any(names(object@hetData) != markers))
	{
		errors <- c(errors, "Slot hetData refers to different markers to slot finals")
	}
	
	#Check that each founder and final is in the pedigree
	if(!all(rownames(object@founders) %in% object@pedigree@lineNames))
	{
		errors <- c(errors, "Not all founder lines were named in the pedigree")
	}
	if(!all(rownames(object@finals) %in% object@pedigree@lineNames))
	{
		errors <- c(errors, "Not all final lines were named in the pedigree")
	}

	#If we have information on the founders in the pedigree, check that the founders ARE in fact those lines
	if(inherits(object@pedigree, "detailedPedigree"))
	{
		if(!all(rownames(object@founders) %in% object@pedigree@lineNames[object@pedigree@initial]) || nrow(object@founders) != length(object@pedigree@initial))
		{
			errors <- c(errors, "Founder lines did not match those specified in the pedigree")
		}
		if(!all(rownames(object@finals) %in% object@pedigree@lineNames[object@pedigree@observed]) || nrow(object@finals) != sum(object@pedigree@observed))
		{
			errors <- c(errors, "Final lines did not match those specified in the pedigree")
		}
	}
	#This checks the relation between the het data, founder data and final data. It doesn't check that the het data is itself valid
	#It also checks that if any of the founders are NULL for a marker, ALL the founder alleles must be NULL, and all the finals alleles must be NULL too
	alleleDataErrors <- .Call("alleleDataErrors", object, 10, PACKAGE="mpMap2")
	if(length(alleleDataErrors) > 0)
	{
		errors <- c(errors, alleleDataErrors)
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.geneticData <- setClass("geneticData", slots=list(finals = "matrix", founders = "matrix", hetData = "hetData", pedigree = "pedigree"), validity = checkGeneticData)
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
					errors <<- c(errors, paste0("Error in geneticData[[", i, "]]: ", geneticErrors))
				}
			}
		}
		errors <- c(errors, checkCompatibleGeneticData(object@geneticData))
	}
	if(!is.null(object@map))
	{
		nMapMarkers <- sum(unlist(lapply(object@map, length)))
		if(nMarkers != nMapMarkers)
		{
			errors <- c(errors, "Number of markers in map is different from the number of markers in slot founders")
		}
		else
		{
			markerNames <- unlist(lapply(object@map, names))
			if(any(colnames(object@founders) != markerNames))
			{
				errors <- c(errors, "Marker names in map disagree with marker names in slot founders")
			}
		}
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.mpcross <- setClass("mpcross", slots = list(geneticData = "list", map = "mapOrNULL"), validity=checkMpcross)

#.compressedMatrix <- setClass("compressedMatrix", slots = list(levels = "numeric", data = "integer"), validity = function(object) .Call("checkCompressedMatrix", object, package="mpMap2"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

checkRF <- function(object)
{
	errors <- c()
	if(storage.mode(object@theta) != "double")
	{
		errors <- c(errors, "storage.mode of slot theta must be double")
	}
	if(!is.null(object@lod) && storage.mode(object@lod) != "double")
	{
		errors <- c(errors, "storage.mode of slot lod must be double")
	}
	if(!is.null(object@lkhd) && storage.mode(object@lkhd) != "double")
	{
		errors <- c(errors, "storage.mode of slot lkhd must be double")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.rf <- setClass("rf", slots = list(r = "numeric", theta = "matrix", lod = "matrixOrNULL", lkhd = "matrixOrNULL"), validity = checkRF)
.mpcrossRF <- setClass("mpcrossRF", contains = "mpcross", slots = list(rf = "rf"))

checkLG <- function(object)
{
	if(!all(object@groups %in% object@allGroups))
	{
		return("Only values in slot allGroups are allowed as values in slot groups")
	}
	return(TRUE)
}
.lg <- setClass("lg", slots = list(groups = "integer", allGroups = "integer"), validity = checkLG)

checkMpcrossLG <- function(object)
{
	if(!is.null(object@map)) return("An mpcross object with assigned linkage groups cannot have a map")
	if(any(names(lg@groups) != colnames(object@founders)))
	{
		return("Marker names implied by names of slots lg@groups and founders were different")
	}
	return(TRUE)
}
.mpcrossLG <- setClass("mpcrossLG", contains = "mpcrossRF", slots = list(lg = "lg"), validity=checkMpcrossLG)