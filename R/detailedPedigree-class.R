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
#' Pedigree for simulation
#' 
#' Class detailedPedigree is similar to the S4 class pedigree, except it also contains information about which lines are going to observed. This allows simulation of a data set with the given pedigree. 
#' @slot initial The indices of the inbred founder lines in the pedigree. These founders lines must be the first lines in the pedigree. 
#' @slot observed A logical vector with one value per line in the pedigree. A value of \code{TRUE} indicates that this line will be genotyped. 
#' @seealso \code{\link[mpMap2]{pedigree-class}}, \code{\link[mpMap2]{simulateMPCross}}, \code{\link[mpMap2]{detailedPedigree}}
#' @include pedigree-class.R
#' @rdname detailedPedigree-class
#' @name  detailedPedigree-class
#' @examples lineNames <- paste0("L", 1:10)
#' mother <- c(0, 0, 1, rep(3, 7))
#' father <- c(0, 0, 2, rep(2, 7))
#' initial <- 1:2
#' lineNames <- paste0("L", 1:10)
#' observed <- c(rep(FALSE, 3), rep(TRUE, 7))
#' detailedPedigreeObj <- detailedPedigree(mother = mother, father = father, initial = initial, 
#' 	observed = observed, lineNames = lineNames, selfing = "finite")
NULL
.detailedPedigree <- setClass("detailedPedigree", contains = "pedigree", slots = list(initial = "integer", observed = "logical"), validity = checkDetailedPedigree)
#' @export
#' @describeIn detailedPedigree-class Construct object of class detailedPedigree
#' @param lineNames The names assigned to the lines.
#' @param mother The female parent of this line, given by name or by index within \code{lineNames}.
#' @param father The male parent of this line, given by name or by index within \code{lineNames}.
#' @param initial The founder lines, given by name or by index within lineNames.
#' @param observed The lines which are observed in the final population, given by name or by index within lineNames.
#' @param selfing Value determining whether or not subsequent analysis of populations generated from this pedigree should assume infinite generations of selfing. Possible values are \code{"finite"} and \code{"infinite"}.
#' @return An object of class \code{detailedPedigree}, suitable for simulation.
#' @seealso \code{\link[mpMap2]{detailedPedigree-class}}
detailedPedigree <- function(lineNames, mother, father, initial, observed, selfing)
{
	lineNames <- as.character(lineNames)
	if(is.character(mother))
	{
		mother <- match(mother, lineNames)
		mother[is.na(mother)] <- 0
	}
	mother <- as.integer(mother)
	if(is.character(father))
	{
		father <- match(father, lineNames)
		father[is.na(father)] <- 0
	}
	father <- as.integer(father)
	if(is.character(initial))
	{
		initial <- match(initial, lineNames)
	}
	initial <- as.integer(initial)
	if(is.character(observed))
	{
		observed <- lineNames %in% observed
	}
	else if(is.numeric(observed))
	{
		observed <- (1:length(lineNames)) %in% observed
	}
	return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = initial, observed = observed, selfing = selfing, warnImproperFunnels = TRUE))
}
