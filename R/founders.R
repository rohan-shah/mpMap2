#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "founders", def = function(object){standardGeneric("founders")})
setMethod(f = "founders", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(object@geneticData[[1]]@founders)
	}
	return(lapply(object@geneticData, function(x) x@founders))
})
setMethod(f = "founders", signature = "geneticData", definition = function(object)
{
	return(object@founders)
})

#' @export
setGeneric(name = "founderNames", def = function(object){standardGeneric("founderNames")})
setMethod(f = "founderNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(rownames(object@geneticData[[1]]@founders))
	}
	return(lapply(object@geneticData, function(x) rownames(x@founders)))
})
setMethod(f = "founderNames", signature = "geneticData", definition = function(object)
{
	return(rownames(object@founders))
})

