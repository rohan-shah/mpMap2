#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Number of genotyped markers
#'
#' Return the number of genotyped markers in an object.
#'
#' If an \code{mpcross} object contains multiple experiments, all experiments are required to have the same markers. So only one number is returned, in all cases. 
#' @rdname nMarkers
#' @param object The \code{mpcross} object from which to extract the marker names
#' @return The number of markers in an object of class \code{mpcross}. 
#' @export
setGeneric(name = "nMarkers", def = function(object){standardGeneric("nMarkers")})
#' @rdname nMarkers
setMethod(f = "nMarkers", signature = "mpcross", definition = function(object)
{
	return(ncol(object@geneticData[[1]]@finals))
})
#' @rdname nMarkers
setMethod(f = "nMarkers", signature = "geneticData", definition = function(object)
{
	return(ncol(object@finals))
})
