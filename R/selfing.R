#' @include pedigree-class.R
#' @include detailedPedigree-class.R
#' @title Get or set a pedigree to have finite or infinite generations of selfing
#' @rdname selfing
#' @description Get or set a pedigree to have finite or infinite generations of selfing
#' @details A pedigree object contains details about the genetic relationships between individuals in a population. Many experiments will include a finite number of generations of inbreeding by selfing, and this information will also be contained in the pedigree. However, when it comes time to actually analyse the poulation, it can be sensible to assume that an infinite number of generations of selfing have actually been performed, as this is computationally quicker. 
#' 
#' This extra information about whether to assume infinite generations of selfing, or the finite number of generations given in the pedigree, is stored in an extra slot, which must have value \code{"finite"} or \code{"infinite"}. If \code{"finite"} is specified, then in subsequent analysis (e.g. computation of IBD genotypes or probabilities) the number of generations of selfing for each line is taken from the pedigree. 
#' @param object The pedigree object for which to get or set the generations of selfing, as finite or infinite. 
#' @param value The new value
#' @return Dimensions of selfing, either \code{"finite"} or \code{"infinite"}.
#' @examples
#' pedigree <- eightParentPedigreeImproperFunnels(initialPopulationSize = 10,
#'      selfingGenerations = 0, nSeeds = 1)
#' selfing(pedigree)
#' selfing(pedigree) <- "finite"
#' @export
setGeneric("selfing<-", function(object, value) standardGeneric("selfing<-"))
#' @export
#' @rdname selfing
setGeneric("selfing", function(object) standardGeneric("selfing"))
#' @rdname selfing
setMethod(f = "selfing", signature = "pedigree", definition = function(object)
{
	object@selfing
})
#' @rdname selfing
setReplaceMethod("selfing", "detailedPedigree", function(object, value)
{
	if(value != "infinite" && value != "finite")
	{
		stop("Selfing must be either \"finite\" or \"infinite\"")
	}
	object@selfing <- value
	object
})
#' @rdname selfing
setReplaceMethod("selfing", "pedigree", function(object, value)
{
	if(value != "infinite" && value != "finite")
	{
		stop("Selfing must be either \"finite\" or \"infinite\"")
	}
	object@selfing <- value
	object
})
