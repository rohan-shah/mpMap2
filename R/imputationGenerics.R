#' @export
#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "flatImputationMapNames", def = function(object){standardGeneric("flatImputationMapNames")})
setMethod(f = "flatImputationMapNames", signature = "imputed", definition = function(x)
{
	return(unlist(lapply(x@map, names)))
})
setMethod(f = "flatImputationMapNames", signature = "geneticData", definition = function(x)
{
	if(is.null(x@imputed))
	{
		stop("Cannot extract the flattened imputation map names if there is no imputation data")
	}
	return(flatImputationMapNames(x@imputed))
})
setMethod(f = "flatImputationMapNames", signature = "mpcrossMapped", definition = function(x)
{
	if(length(x@geneticData) == 1)
	{
		return(flatImputationMapNames(x@geneticData[[1]]))
	}
	else
	{
		return(lapply(x@geneticData, flatImputationMapNames))
	}
})
#' @export
setGeneric(name = "imputationMap", def = function(object){standardGeneric("imputationMap")})
setMethod(f = "imputationMap", signature = "imputed", definition = function(x)
{
	return(x@map)
})
setMethod(f = "imputationMap", signature = "geneticData", definition = function(x)
{
	if(is.null(x@imputed))
	{
		stop("Cannot extract the imputation map if there is no imputation data")
	}
	return(imputationMap(x@imputed))
})
setMethod(f = "imputationMap", signature = "mpcrossMapped", definition = function(x)
{
	if(length(x@geneticData) == 1)
	{
		return(imputationMap(x@geneticData[[1]]))
	}
	else
	{
		return(lapply(x@geneticData, imputationMap))
	}
})
#' @export
setGeneric(name = "imputationData", def = function(object){standardGeneric("imputationData")})
setMethod(f = "imputationData", signature = "imputed", definition = function(x)
{
	return(x@data)
})
setMethod(f = "imputationData", signature = "geneticData", definition = function(x)
{
	if(is.null(x@imputed))
	{
		stop("Imputation data is not present")
	}
	return(imputationData(x@imputed))
})
setMethod(f = "imputationData", signature = "mpcrossMapped", definition = function(x)
{
	if(length(x@geneticData) == 1)
	{
		return(imputationData(x@geneticData[[1]]))
	}
	else
	{
		return(lapply(x@geneticData, imputationData))
	}
})
