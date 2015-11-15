checkPedigree <- function(object)
{
	nTotalLines <- length(object@lineNames)
	errors <- c()
	if(length(object@mother) != nTotalLines || length(object@father) != nTotalLines)
	{
		errors <- c(errors, "Lengths of slots lineNames, mother and father must be the same")
	}
	if(length(errors) > 0) return(errors)
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

	#Line names must be unique
	if(length(unique(object@lineNames)) != nTotalLines)
	{
		errors <- c(errors, "Line names must be unique")
	}
	if(nFounders(object) < 2)
	{
		errors <- c(errors, "Number of founders must be at least 2")
	}

	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.pedigree <- setClass("pedigree", slots = list(lineNames = "character", mother = "integer", father = "integer", selfing = "character"), validity = checkPedigree)
