#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Genetic data for founding lines
#'
#' Return the genetic data matrix for the founding lines
#'
#' If the \code{mpcross} object contains a single experiment a matrix is returned, with rows corresponding to founding lines and columns corresponding to markers. If an \code{mpcross} object contains multiple experiments a list of such matrices is returned, one for each experiment.
#' @rdname founders
#' @param object The \code{mpcross} object from which to extract the genetic data matrix of the founding lines
#' @return An integer matrix, with rows corresponding to founding lines and columns corresponding to markers, or a list of such matrices in the case of multiple experiments. 
#' @export
setGeneric(name = "founders", def = function(object){standardGeneric("founders")})
#' @rdname founders
setMethod(f = "founders", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(object@geneticData[[1]]@founders)
	}
	return(lapply(object@geneticData, function(x) x@founders))
})
#' @rdname founders
setMethod(f = "founders", signature = "geneticData", definition = function(object)
{
	return(object@founders)
})
#' @title Names of founding genetic lines
#'
#' Return the names of the founding genetic lines
#'
#' If the \code{mpcross} object contains a single experiment a vector of names of genetic lines is returned. If an \code{mpcross} object contains multiple experiments a list of vectors of names is returned. 
#' @param object The \code{mpcross} object from which to extract the names of the founding genetic lines
#' @return A vector of names of genetic lines, or a list of such vectors, in the case of multiple experiments. 
#' @rdname founderNames
#' @export
setGeneric(name = "founderNames", def = function(object){standardGeneric("founderNames")})
#' @rdname founderNames
setMethod(f = "founderNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(rownames(object@geneticData[[1]]@founders))
	}
	return(lapply(object@geneticData, function(x) rownames(x@founders)))
})
#' @rdname founderNames
setMethod(f = "founderNames", signature = "geneticData", definition = function(object)
{
	return(rownames(object@founders))
})

