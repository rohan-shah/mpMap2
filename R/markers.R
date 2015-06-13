setGeneric(name = "markers", def = function(object){standardGeneric("markers")})
setMethod(f = "markers", signature = "mpcross", definition = function(object)
{
	return(markers(object@geneticData[[1]]))
})
setMethod(f = "markers", signature = "geneticData", definition = function(object)
{
	return(colnames(object@finals))
})
setMethod(f = "markers", signature = "rf", definition = function(object)
{
	return(colnames(object@theta))
})
setMethod(f = "markers", signature = "lg", definition = function(object)
{
	return(names(object@groups))
})
setMethod(f = "markers", signature = "hetData", definition = function(object)
{
	return(names(object))
})