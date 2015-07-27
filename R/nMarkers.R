#' @include mpcross-class.R
#' @include geneticData-class.R
setGeneric(name = "nMarkers", def = function(object){standardGeneric("nMarkers")})
setMethod(f = "nMarkers", signature = "mpcross", definition = function(object)
{
	return(ncol(object@geneticData[[1]]@finals))
})
setMethod(f = "nMarkers", signature = "geneticData", definition = function(object)
{
	return(ncol(object@finals))
})