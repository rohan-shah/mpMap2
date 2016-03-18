#' Pedigree class
#' 
#' This class describes a pedigree for an experimental design. Although package mpMap2 only allows for the analysis of pedigrees corresponding to multi-parent crosses, this pedigree class can describe arbitrary experimental designs. 
#' @slot mother The index within the pedigree of the mother of this individual
#' @slot father The index within the pedigree of the father of this individual
#' @slot lineName The name of this individual
#' @slot selfing A value indicating whether analysis of an experiment using this pedigree should assume infinite generations of selfing. A value of \code{"infinite"} indicates infinite generations of selfing, and a value of \code{"finite"} indicates finite generations of selfing. 
#' @seealso \code{\link[mpMap2]{pedigree-class}}, \code{\link[mpMap2]{simulateMPCross}}, \code{\link[mpMap2]{detailedPedigree-class}}, \code{\link[mpMap2]{rilPedigree}}, \code{\link[mpMap2]{f2Pedigree}}, \code{\link[mpMap2]{fourParentPedigreeRandomFunnels}}, \code{\link[mpMap2]{fourParentPedigreeSingleFunnel}}, \code{\link[mpMap2]{eightParentPedigreeRandomFunnels}}, \code{\link[mpMap2]{eightParentPedigreeSingleFunnel}}, \code{\link[mpMap2]{sixteenParentPedigreeRandomFunnels}}
#' @include pedigree-class.R
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
	#Entry selfing must be either "finite" of "infinite"
	if(object@selfing != "finite" && object@selfing != "infinite")
	{
		errors <- c(errors, "Slot selfing must be either \"infinite\" or \"finite\"")
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
