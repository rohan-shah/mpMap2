setGeneric(name = "markers", def = function(object){standardGeneric("markers")})
setMethod(f = "markers", signature = "mpcross", definition = function(object)
{
	return(markers(object@geneticData[[1]]))
})
setMethod(f = "markers", signature = "geneticData", definition = function(object)
{
	return(colnames(object@finals))
})