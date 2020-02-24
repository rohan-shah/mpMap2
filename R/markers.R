#' @include mpcross-class.R
#' @include rf-class.R
#' @include lg-class.R
#' @include hetData-class.R
#' @title Genotyped markers
#'
#' Return the names of the genotyped markers. 
#'
#' If an \code{mpcross} object contains multiple experiments, all experiments are required to have the same markers. So a single vector of marker names is returned, in all cases. 
#' @rdname markers
#' @param object The \code{mpcross} object from which to extract the marker names
#' @return The names of the genetic markers.
#' @export
setGeneric(name = "markers", def = function(object){standardGeneric("markers")})
#' @rdname markers
setMethod(f = "markers", signature = "mpcross", definition = function(object)
{
	return(markers(object@geneticData[[1]]))
})
#' @rdname markers
setMethod(f = "markers", signature = "geneticData", definition = function(object)
{
	return(colnames(object@finals))
})
#' @rdname markers
setMethod(f = "markers", signature = "rf", definition = function(object)
{
	return(object@theta@markers)
})
#' @rdname markers
setMethod(f = "markers", signature = "lg", definition = function(object)
{
	return(names(object@groups))
})
#' @rdname markers
setMethod(f = "markers", signature = "hetData", definition = function(object)
{
	return(names(object))
})
