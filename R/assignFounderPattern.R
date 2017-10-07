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
	if(ncol(e2@data) != nMarkers(e1))
	{
		stop("Input founder matrix had the wrong number of markers")
	}
	if(any(sort(markers(e1)) != sort(colnames(e2@data))))
	{
		stop("Input founder matrix had different markers to the dataset")
	}
	if(any(markers(e1) != colnames(e2@data)))
	{
		warning("Markers in the founder matrix were in a different order to the markers in the dataset; reordering the columns of the founder matrix")
		e2@data <- e2@data[,markers(e1)]
	}
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
	if(ncol(e2@data) != nMarkers(e1))
	{
		stop("Input founder matrix had the wrong number of markers")
	}
	if(any(sort(markers(e1)) != sort(colnames(e2@data))))
	{
		stop("Input founder matrix had different markers to the dataset")
	}
	if(any(markers(e1) != colnames(e2@data)))
	{
		warning("Markers in the founder matrix were in a different order to the markers in the dataset; reordering the columns of the founder matrix")
		e2@data <- e2@data[,markers(e1)]
	}
	if(length(e1@geneticData) != 1)
	{
		stop("Can only apply assignFounderPattern to an mpcross object with a single experiment")
	}
	result <- .Call("assignFounderPattern", e1@geneticData[[1]], founderPattern = e2@data)
	geneticData <- new("geneticData", founders = result@founders, finals = result@finals, pedigree = e1@geneticData[[1]]@pedigree, map = e1@map, hetData = result@hetData)
	geneticDataList <- new("geneticDataList", list(geneticData))
	return(new("mpcrossMapped", geneticData = geneticDataList,  map = e1@map))
})
