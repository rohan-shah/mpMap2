#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
setGeneric(name = "flatImputationMapNames", def = function(object, ...){standardGeneric("flatImputationMapNames")})
setMethod(f = "flatImputationMapNames", signature = "imputed", definition = function(object, ...)
{
	return(unlist(lapply(object@map, names)))
})
setMethod(f = "flatImputationMapNames", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Cannot extract the flattened imputation map names if there is no imputation data")
	}
	return(flatImputationMapNames(object@imputed))
})
setMethod(f = "flatImputationMapNames", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("design" %in% names(args) && args$design != 1) stop("Input design specified and not equal to one, but there is only one design in this object")
		return(flatImputationMapNames(object@geneticData[[1]]))
	}
	else
	{
		if("design" %in% names(args))
		{
			return(flatImputationMapNames(object@geneticData[[args$design]]))
		}
		else
		{
			return(lapply(object@geneticData, flatImputationMapNames))
		}
	}
})
#' @export
setGeneric(name = "imputationMap", def = function(object, ...){standardGeneric("imputationMap")})
setMethod(f = "imputationMap", signature = "imputed", definition = function(object, ...)
{
	return(object@map)
})
setMethod(f = "imputationMap", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Cannot extract the imputation map if there is no imputation data")
	}
	return(imputationMap(object@imputed))
})
setMethod(f = "imputationMap", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("design" %in% names(args) && args$design != 1) stop("Input design specified and not equal to one, but there is only one design in this object")
		return(imputationMap(object@geneticData[[1]]))
	}
	else
	{
		if("design" %in% names(args))
		{
			return(imputationMap(object@geneticData[[args$design]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationMap))
		}
	}
})
#' @export
setGeneric(name = "imputationData", def = function(object, ...){standardGeneric("imputationData")})
setMethod(f = "imputationData", signature = "imputed", definition = function(object, ...)
{
	return(object@data)
})
setMethod(f = "imputationData", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Imputation data is not present")
	}
	return(imputationData(object@imputed))
})
setMethod(f = "imputationData", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("design" %in% names(args) && args$design != 1) stop("Input design specified and not equal to one, but there is only one design in this object")
		return(imputationData(object@geneticData[[1]]))
	}
	else
	{
		if("design" %in% names(args))
		{
			return(imputationData(object@geneticData[[args$design]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationData))
		}
	}
})
#' @export
setGeneric(name = "imputationKey", def = function(object, ...){standardGeneric("imputationKey")})
setMethod(f = "imputationKey", signature = "imputed", definition = function(object, ...)
{
	return(object@key)
})
setMethod(f = "imputationKey", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Imputation data is not present")
	}
	return(imputationKey(object@imputed))
})
setMethod(f = "imputationKey", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("design" %in% names(args) && args$design != 1) stop("Input design specified and not equal to one, but there is only one design in this object")
		return(imputationKey(object@geneticData[[1]]))
	}
	else
	{
		if("design" %in% names(args))
		{
			return(imputationKey(object@geneticData[[args$design]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationKey))
		}
	}
})
