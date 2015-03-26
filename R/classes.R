checkPedigree <- function(object)
{
	nLines <- length(object@lineNames)
	errors <- c()
	if(length(object@mother) != nLines || length(object@father) != nLines || length(object@observed) != nLines)
	{
		errors <- c(errors, "Lengths of slots lineNames, mother, father and observed must be the same")
	}

	if(any(is.na(object@initial)))
	{
		errors <- c("errors", "Slot initial contained NA values")
	}
	if(any(is.na(object@lineNames)))
	{
		errors <- c("errors", "Slot lineNames contained NA values")
	}
	if(any(is.na(object@mother)))
	{
		errors <- c("errors", "Slot mother contained NA values")
	}
	if(any(is.na(object@father)))
	{
		errors <- c("errors", "Slot father contained NA values")
	}
	if(any(is.na(object@observed)))
	{
		errors <- c("errors", "Slot observed contained NA values")
	}

	if(length(object@initial) == 0)
	{
		errors <- c(errors, "Slot initial must be non-empty")
	}

	if(any(object@initial < 0 | object@initial > nLines))
	{
		errors <- c(errors, "Values in slot initial had invalid values")
	}
	if(any(object@mother < 0 | object@mother > nLines))
	{
		errors <- c(errors, "Values in slot mother had invalid values")
	}
	if(any(object@father < 0 | object@father > nLines))
	{
		errors <- c(errors, "Values in slot father had invalid values")
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

	#Parents must come first in the pedigree
	if(!all(object@mother < 1:nLines) || !all(object@father < 1:nLines))
	{
		errors <- c(errors, "Mother and father must preceed offspring in the pedigree")
	}

	if(!all(object@initial == 1:max(object@initial)))
	{
		errors <- c(errors, "Initial lines must be at the start of the pedigree")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
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
.pedigree <- setClass("pedigree", slots = list(lineNames = "character", mother = "integer", father = "integer", initial = "integer", observed = "logical"), validity = checkPedigree)
checkMpcross <- function(object)
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

	nMarkers <- ncol(object@founders)
	if(ncol(object@finals) != nMarkers || length(object@hetData) != nMarkers)
	{
		errors <- c(errors, "Slots finals, founders and hetData had different numbers of markers")
	}
	if(any(colnames(object@founders) != colnames(object@finals)))
	{
		errors <- c(errors, "Slot finals must have the same colnames as slot founders")
	}
	if(any(names(object@hetData) != colnames(object@finals)))
	{
		errors <- c(errors, "Slot hetData refers to different markers to slot finals")
	}

	if(nrow(object@finals) != sum(pedigree@observed))
	{
		errors <- c(errors, "Number of rows of slot finals was inconsistent with slot observed of pedigree")
	}
	if(any(rownames(object@finals) != pedigree@lineNames[pedigree@observed]))
	{
		errors <- c(errors, "Row names of slot finals were inconsistent with slot observed of pedigree")
	}
	alleleDataErrors <- .Call("alleleDataErrors", object, PACKAGE="mpMap2")
	if(length(alleleDataErrors) > 0)
	{
		if(length(alleleDataErrors) > 10)
		{
			errors <- c(alleleDataErrors, "Not reporting futher errors related to invalid alleles", errors)
		}
		else errors <- c(alleleDataErrors, errors)
	}
	return(errors)
}
.mpcross <- setClass("mpcross", slots = list(finals = "matrix", founders = "matrix", hetData = "hetData", pedigree = "pedigree"), validity=checkMpcross)