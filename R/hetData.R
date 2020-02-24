#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Get the encoding of marker heterozygotes
#' @description Get the encoding of marker heterozygotes
#' @rdname hetData
#' @details Get the encoding of markers heterozygotes, either for all markers, or a specific marker. 
#' @param object The object from which to extract the encoding data
#' @param marker The marker of interest. If this is missing, heterozygote encoding data is returned for all markers. 
#' @return Heterozygote encoding data, for either a specific marker or all markers. 
#' @export
setGeneric(name = "hetData", def = function(object, marker){standardGeneric("hetData")})
#' @rdname hetData
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
#' @rdname hetData
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
