#' @include detailedPedigree-class.R
#' @include pedigree-class.R
#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Number of genotyped markers
#'
#' Return the number of genotyped markers in an object.
#'
#' If an \code{mpcross} object contains multiple experiments, one number is returned per experiment. 
#' @rdname nFounders
#' @param object The \code{mpcross} object from which to extract the number of founders
#' @return The number of founding lines in the population, or a list of numbers in the case of multiple experiments contained in a single object.  
#' @export
setGeneric(name = "nFounders", def = function(object){standardGeneric("nFounders")})
#' @rdname nFounders
setMethod(f = "nFounders", signature = "detailedPedigree", definition = function(object)
{
	return(length(object@initial))
})
#' @rdname nFounders
setMethod(f = "nFounders", signature = "pedigree", definition = function(object)
{
	return(sum(object@mother == 0 & object@father == 0))
})
#' @rdname nFounders
setMethod(f = "nFounders", signature = "mpcross", definition = function(object)
{
	return(unlist(lapply(object@geneticData, nFounders)))
})
#' @rdname nFounders
setMethod(f = "nFounders", signature = "geneticData", definition = function(object)
{
	return(nrow(object@founders))
})
