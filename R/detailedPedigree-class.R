#' Pedigree for simulation
#' 
#' Class detailedPedigree is similar to the S4 class pedigree, except it also contains information about which lines are going to observed. This allows us to simulate a data set with the given pedigree. 
#' @slot initial The indices of the inbred founder lines in the pedigree. These founders lines must be the first lines in the pedigree. 
#' @slot observed A logical vector with one value per line in the pedigree. A value of \code{TRUE} indicates that this line will be genotyped. 
#' @seealso \code{\link[mpMap2]{pedigree-class}}, \code{\link[mpMap2]{simulateMPCross}}
#' @include pedigree-class.R
#' @example
#' 
checkDetailedPedigree <- function(object)
{
	nTotalLines <- length(object@lineNames)
	errors <- c()
	#length(object@mother) == length(object@lineNames) and length(object@father) == length(object@lineNames) are checked in checkPedigree
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

	if(any(object@initial <= 0 | object@initial > nTotalLines))
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
	if(length(errors) > 0) return(errors)

	if(!all(sort(object@initial) == 1:max(object@initial)))
	{
		errors <- c(errors, "Initial lines must be at the start of the pedigree")
	}
	#Check that everything with mother or father zero is in the initial slot
	if(any(object@mother[max((object@initial)+1):length(object@mother)] == 0) || any(object@father[(max(object@initial)+1):length(object@father)] == 0))
	{
		errors <- c(errors, "Some founder lines were not listed in slot initial")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.detailedPedigree <- setClass("detailedPedigree", contains = "pedigree", slots = list(initial = "integer", observed = "logical"), validity = checkDetailedPedigree)