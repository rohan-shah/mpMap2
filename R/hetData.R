#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "hetData", def = function(object, marker){standardGeneric("hetData")})
setMethod(f = "hetData", signature = c("mpcross", "ANY"), definition = function(object, marker)
{
	if(length(object@geneticData) == 1)
	{
		if(missing(marker))
		{
			return(object@geneticData[[1]]@hetData)
		}
		else 
		{
			return(object@geneticData[[1]]@hetData[[marker]])
		}
	}
	return(lapply(object@geneticData, function(x) x@hetData[[marker]]))
})
setMethod(f = "hetData", signature = c("geneticData", "ANY"), definition = function(object, marker)
{
	if(missing(marker))
	{
		return(object@hetData)
	}
	else
	{
		return(object@hetData[[marker]])
	}
})
