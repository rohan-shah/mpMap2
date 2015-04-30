setGeneric(name = "nFounders", def = function(object){standardGeneric("nFounders")})
setMethod(f = "nFounders", signature = "pedigree", definition = function(object)
{
	return(length(object@initial))
})
setMethod(f = "nFounders", signature = "mpcross", definition = function(object)
{
	return(unlist(lapply(object@geneticData, nFounders)))
})
setMethod(f = "nFounders", signature = "geneticData", definition = function(object)
{
	return(nrow(object@founders))
})