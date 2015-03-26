setGeneric(name = "nFounders", def = function(object){standardGeneric("nFounders")})
setMethod(f = "nFounders", signature = "pedigree", definition = function(object)
{
	return(length(object@initial))
})
setMethod(f = "nFounders", signature = "mpcross", definition = function(object)
{
	return(length(object@pedigree@initial))
})