#' @include detailedPedigree-class.R
#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Number of genotyped lines
#'
#' Return the number of genotyped lines in an object.
#'
#' This includes only the number of final lines genotyped in the population, and does not include the founding lines. If an \code{mpcross} object contains multiple experiments, one number is returned per experiment. 
#' @param object The \code{mpcross} object from which to extract the number of genotyped lines.
#' @return The number of genetic lines in the population, or a list of numbers in the case of multiple experiments contained in a single object.  
#' @rdname nLines
#' @export
setGeneric(name = "nLines", def = function(object){standardGeneric("nLines")})
#' @rdname nLines
setMethod(f = "nLines", signature = "mpcross", definition = function(object)
{
	return(unlist(lapply(object@geneticData, nLines)))
})
#' @rdname nLines
setMethod(f = "nLines", signature = "geneticData", definition = function(object)
{
	return(nrow(object@finals))
})
