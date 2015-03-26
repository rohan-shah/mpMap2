setGeneric(name = "nLines", def = function(object){standardGeneric("nLines")})
setMethod(f = "nLines", signature = "pedigree", definition = function(object)
{
	return(length(object@initial))
})
setMethod(f = "nLines", signature = "mpcross", definition = function(object)
{
	return(nrow(object@finals))
})