#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Get IBD probability data
#' @description Get the identity-by-descent probability data from an mpcross object.
#' @details mpMap2 stores IBD probabilities in a matrix where the number of rows is the number of alleles times the number of genetic lines, and the columns are the positions at which probabilities are calculated. In the example below, the row names are \code{L115 - L1, L115 - L2, L115 - L3, L115 - L4, L120 - L1}, etc, and the column names are \code{DXM1, DXM2, DXM3}, etc. So, for example, for a population generated from four founders and assumed to be totally inbred, the first four values in the first column are the probabilities that genetic line 1 carries alleles from specific founders, at a specific position. The first four columns give the probabilities for genetic line 2 at the next position, etc. 
#'
#' This can be an inconvenient layout for some operations. This function returns a matrix where the alleles appear as part of the columns, rather than the rows. For example, after applying this function to the given example, the first four values in the first row will be the probabilities that genetic line 1 carries alleles from specific founders, at a specific position. 
#' @rdname probabilityData
#' @param object The \code{mpcross} object from which to extract the probability data.
#' @param ... Additional options. Only \code{design} is supported, and gets the probability data for only a single experiment. 
#' @return A numeric matrix containing the IBD probabliity data, or a list of such matrices in the case of multiple experiments within a single object. 
#' @examples
#' data(simulatedFourParentData)
#' crossSNP <- simulatedFourParentData + multiparentSNP(keepHets = FALSE)
#' mapped <- mpcrossMapped(crossSNP, map = simulatedFourParentMap)
#' probabilities <- computeGenotypeProbabilities(mapped, error = 0.05)
#' probabilityData <- probabilityData(probabilities)
#' probabilityData[1:5, 1:5]
#' @export
setGeneric(name = "probabilityData", def = function(object, ...){standardGeneric("probabilityData")})
#' @rdname probabilityData
setMethod(f = "probabilityData", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@probabilities))
	{
		stop("Cannot extract the probabilities if there is no probability data")
	}
	return(object@probabilities@data)
})
#' @rdname probabilityData
setMethod(f = "probabilityData", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("design" %in% names(args) && args$design != 1) stop("Input design specified and not equal to one, but there is only one design in this object")
		return(probabilityData(object@geneticData[[1]]))
	}
	else
	{
		if("design" %in% names(args))
		{
			return(probabilityData(object@geneticData[[args$design]]))
		}
		else
		{
			return(lapply(object@geneticData, probabilityData))
		}
	}
})

