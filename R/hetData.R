#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "hetData", def = function(object){standardGeneric("hetData")})
setMethod(f = "hetData", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(object@geneticData[[1]]@hetData)
	}
	return(lapply(object@geneticData, function(x) x@hetData))
})
setMethod(f = "hetData", signature = "geneticData", definition = function(object)
{
	return(object@hetData)
})