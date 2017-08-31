#' @include mpcross-class.R
#' @include geneticData-class.R
setClass("assignFounderPattern", slots = list(data = "matrix"))
#' @export
assignFounderPattern <- function(founderMatrix)
{
	return(new("assignFounderPattern", data = founderMatrix))
}
#' @rdname internalOperators
#' @title Internal operators for mpMap2
setMethod(f = "+", signature = c("geneticData", "assignFounderPattern"), definition = function(e1, e2)
{
	result <- .Call("assignFounderPattern", e1, founderPattern = e2@data)
	return(result)
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "assignFounderPattern"), definition = function(e1, e2)
{
	if(length(e1@geneticData) != 1)
	{
		stop("Can only apply assignFounderPattern to an mpcross object with a single experiment")
	}
	result <- .Call("assignFounderPattern", e1@geneticData[[1]], founderPattern = e2@data)
	geneticData <- new("geneticData", founders = result@founders, finals = result@finals, pedigree = e1@geneticData[[1]]@pedigree, hetData = result@hetData)
	geneticDataList <- new("geneticDataList", list(geneticData))
	return(new("mpcross", geneticData = geneticDataList))
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcrossMapped", "assignFounderPattern"), definition = function(e1, e2)
{
	if(length(e1@geneticData) != 1)
	{
		stop("Can only apply assignFounderPattern to an mpcross object with a single experiment")
	}
	result <- .Call("assignFounderPattern", e1@geneticData[[1]], founderPattern = e2@data)
	geneticData <- new("geneticData", founders = result@founders, finals = result@finals, pedigree = e1@geneticData[[1]]@pedigree, map = e1@map, hetData = result@hetData)
	geneticDataList <- new("geneticDataList", list(geneticData))
	return(new("mpcrossMapped", geneticData = geneticDataList,  map = e1@map))
})
