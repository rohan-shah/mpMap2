setGeneric(name = "nLines", def = function(object){standardGeneric("nLines")})
setMethod(f = "nLines", signature = "detailedPedigree", definition = function(object)
{
	return(sum(object@observed))
})
setMethod(f = "nLines", signature = "mpcross", definition = function(object)
{
	return(unlist(lapply(object@geneticData, nLines)))
})
setMethod(f = "nLines", signature = "geneticData", definition = function(object)
{
	return(nrow(object@finals))
})