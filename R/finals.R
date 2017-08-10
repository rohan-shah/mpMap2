#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "finals", def = function(object){standardGeneric("finals")})
setMethod(f = "finals", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(object@geneticData[[1]]@finals)
	}
	return(lapply(object@geneticData, function(x) x@finals))
})
setMethod(f = "finals", signature = "geneticData", definition = function(object)
{
	return(object@finals)
})

#' @export
setGeneric(name = "finalNames", def = function(object){standardGeneric("finalNames")})
setMethod(f = "finalNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1)
	{
		return(rownames(object@geneticData[[1]]@finals))
	}
	return(lapply(object@geneticData, function(x) rownames(x@finals)))
})
setMethod(f = "finalNames", signature = "geneticData", definition = function(object)
{
	return(rownames(object@finals))
})
