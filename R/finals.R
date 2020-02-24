#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Genetic data for final lines
#'
#' Return the genetic data matrix for the final lines
#'
#' If the \code{mpcross} object contains a single experiment a matrix is returned, with rows corresponding to genotyped lines and columns corresponding to markers. The founding lines of the population are excluded from this matrix. If an \code{mpcross} object contains multiple experiments a list of such matrices is returned, one for each experiment.
#' @rdname finals
#' @param object The \code{mpcross} object from which to extract the genetic data matrix 
#' @return An integer matrix with rows corresponding to genotyped lines and columns corresponding to markers.
#' @export
setGeneric(name = "finals", def = function(object){standardGeneric("finals")})
#' @rdname finals
setMethod(f = "finals", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(object@geneticData[[1]]@finals)
	}
	return(lapply(object@geneticData, function(x) x@finals))
})
#' @rdname finals
setMethod(f = "finals", signature = "geneticData", definition = function(object)
{
	return(object@finals)
})
#' @title Names of genetic lines
#'
#' Return the names of the genetic lines
#'
#' If the \code{mpcross} object contains a single experiment a vector of names of genetic lines is returned. The names of the founding lines for the population are excluded. If an \code{mpcross} object contains multiple experiments a list of vectors of names is returned. 
#' @param object The \code{mpcross} object from which to extract the names of the genetic lines
#' @return The names of the genetic lines in the final population. 
#' @rdname finalNames
#' @export
setGeneric(name = "finalNames", def = function(object){standardGeneric("finalNames")})
#' @rdname finalNames
setMethod(f = "finalNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(rownames(object@geneticData[[1]]@finals))
	}
	return(lapply(object@geneticData, function(x) rownames(x@finals)))
})
#' @rdname finalNames
setMethod(f = "finalNames", signature = "geneticData", definition = function(object)
{
	return(rownames(object@finals))
})
