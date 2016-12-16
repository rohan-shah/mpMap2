#' @include hetData-class.R
#' @include pedigree-class.R
#' @include map-class.R
checkGeneticData <- function(object)
{
	errors <- c()
	if(!is.numeric(object@founders))
	{
		errors <- c(errors, "Slot founders must be a numeric matrix")
	}
	#Error if all founders are NA (because the next condition after this involves max(,na.rm=T), which requires at least one non-NA val)
	else if(all(is.na(object@founders)))
	{
		errors <- c(errors, "Slot founders cannot contain only NA")
	}
	else if(max(abs(round(object@founders) - object@founders), na.rm=TRUE) > 0)
	{
		errors <- c(errors, "Slot founders must contain integer values")
	}
	
	if(!is.numeric(object@finals))
	{
		errors <- c(errors, "Slot finals must be a numeric matrix")
	}
	#Don't do the next two checks if there is no data.
	else if(length(object@finals) > 0)
	{
		if(all(is.na(object@finals)))
		{
			errors <- c(errors, "Slot finals cannot contain only NA")
		}
		else if(max(abs(round(object@finals) - object@finals), na.rm=TRUE) > 0)
		{
			errors <- c(errors, "Slot finals must contain integer values")
		}
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
	if(is.null(dimnames(object@founders)) || any(unlist(lapply(dimnames(object@founders), is.null))))
	{
		errors <- c(errors, "Slot founders must have row and column names")
	}
	#Don't check for row and column names if there is no actual data.
	if(length(object@finals) > 0 && (is.null(dimnames(object@finals)) || any(unlist(lapply(dimnames(object@finals), is.null)))))
	{
		errors <- c(errors, "Slot finals must have row and column names")
	}

	#Check for NA's in row and column names of finals and founders
	if(any(unlist(lapply(dimnames(object@founders), function(x) any(is.na(x))))))
	{
		errors <- c(errors, "Row and column names of slot founders cannot be NA")
	}
	#Skip this if there is no data (this avoids a warnings)
	if(length(object@finals) > 0 && any(unlist(lapply(dimnames(object@finals), function(x) any(is.na(x))))))
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
		#If they have the wronh dimensions then the colnames / names checks have the wrong lengths. So just return from here. 
		return(errors)
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
		#Only check this if there are actually final lines
		if(length(object@finals) > 0)
		{
			if(!all(rownames(object@finals) %in% object@pedigree@lineNames[object@pedigree@observed]) || nrow(object@finals) != sum(object@pedigree@observed))
			{
				errors <- c(errors, "Final lines did not match those specified in the pedigree")
			}
		}
	}
	#Don't continue to the C code if there are already problems - Especially if the founders / finals slots have the wrong types
	if(length(errors) > 0) return(errors)
	#This checks the relation between the het data, founder data and final data. It doesn't check that the het data is itself valid
	#It also checks that if any of the founders are NULL for a marker, ALL the founder alleles must be NULL, and all the finals alleles must be NULL too
	alleleDataErrors <- .Call("alleleDataErrors", object, 10, PACKAGE="mpMap2")
	if(length(alleleDataErrors) > 0)
	{
		errors <- c(errors, alleleDataErrors)
	}
	if(length(errors) > 0) return(errors)

	#Check imputed slot
	if(!is.null(object@imputed))
	{
		if(!identical(dim(object@imputed@data), dim(object@finals)))
		{
			return("Dimensions of slot imputed@data must be the same as those of slot finals")
		}
		if(!identical(dimnames(object@imputed@data), dimnames(object@finals)))
		{
			return("Row and column names of slot imputed@data must be the same as those of slot finals")
		}
		errors <- validObject(object@imputed)
		if(length(errors) > 0) return(errors)
	}
	#Check probabilities
	if(!is.null(object@probabilities))
	{
		if(!is.numeric(object@probabilities@data))
		{
			return("Slot probabilities@data must be a numeric matrix")
		}
		nGenotypes <- length(unique(object@probabilities@key[,3]))
		if(nrow(object@probabilities@data) != nLines(object) * nGenotypes)
		{
			return("Number of rows of probabilities@data must be consistent with probabilities@key and nrow(finals)")
		}
		if(!identical(colnames(object@probabilities@data), markers(object)))
		{
			return("Object probabilities@data had the wrong column names")
		}
		errors <- validObject(object@probabilities)
		if(length(errors) > 0) return(errors)
	}
	return(TRUE)
}
checkImputedData <- function(object)
{
	if(!is.numeric(object@data))
	{
		return("Slot data must be an integer matrix")
	}
	if(!is.numeric(object@key))
	{
		return("Slot key must be an integer matrix")
	}
	if(ncol(object@key) != 3L)
	{
		return("Slot key must have three columns")
	}
	allMapMarkers <- unlist(lapply(object@map, names))
	if(!isTRUE(all.equal(allMapMarkers, colnames(object@data))))
	{
		return("Slot data must have marker names that match the markers in slot map")
	}
	if(!.Call("checkImputedBounds", object, PACKAGE="mpMap2"))
	{
		return("Slot imputed@data must contain values in imputed@key")
	}
}
.imputed <- setClass("imputed", slots=list(data = "matrix", key = "matrix", map = "map"), validity = checkImputedData)
setClassUnion("imputedOrNULL", c("imputed", "NULL"))
checkProbabilities <- function(object)
{
	if(!is.numeric(object@data))
	{
		return("Slot data must be an integer matrix")
	}
	if(!is.numeric(object@key))
	{
		return("Slot key must be an integer matrix")
	}
	if(ncol(object@key) != 3L)
	{
		return("Slot key must have three columns")
	}
}
.probabilities <- setClass("probabilities", slots=list(data = "matrix", key = "matrix", map = "map"), validity = checkProbabilities)
setClassUnion("probabilitiesOrNULL", c("probabilities", "NULL"))
.geneticData <- setClass("geneticData", slots=list(finals = "matrix", founders = "matrix", hetData = "hetData", pedigree = "pedigree", imputed = "imputedOrNULL", probabilities = "probabilitiesOrNULL"), validity = checkGeneticData)
checkGeneticDataList <- function(object)
{
	if(any(unlist(lapply(object, class)) != "geneticData"))
	{
		return("An object of class geneticDataList must contain objects of class geneticData")
	}
	errors <- lapply(object, validObject)
	if(sum(unlist(lapply(errors, is.logical))) == length(object))
	{
		return(TRUE)
	}
	sapply(1:length(object), function(x) errors[[x]] <<- paste0("Error in geneticData object ", x, ": ", errors[[x]]))
	return(unlist(errors))
}
.geneticDataList <- setClass("geneticDataList", "list", validity = checkGeneticDataList)
