#' @title Transpose IBD probabilities
#' @description Transpose the IBD probabilities matrix, so that the different alleles or founders appear on the columns, rather than the rows
#' @details mpMap2 stores IBD probabilities in a matrix where the number of rows is the number of alleles times the number of genetic lines, and the columns are the positions at which probabilities are calculated. In the example below, the row names are \code{L115 - L1, L115 - L2, L115 - L3, L115 - L4, L120 - L1}, etc, and the column names are \code{DXM1, DXM2, DXM3}, etc. So, for example, for a population generated from four founders and assumed to be totally inbred, the first four values in the first column are the probabilities that genetic line 1 carries alleles from specific founders, at a specific position. The first four columns give the probabilities for genetic line 2 at the next position, etc. 
#'
#' This can be an inconvenient layout for some operations. This function returns a matrix where the alleles appear as part of the columns, rather than the rows. For example, after applying this function to the given example, the first four values in the first row will be the probabilities that genetic line 1 carries alleles from specific founders, at a specific position. 
#' @param inputObject The \code{mpcross} object containing the probability data. 
#' @return A numeric matrix containing IBD probability data. 
#' @examples
#' data(simulatedFourParentData)
#' crossSNP <- simulatedFourParentData + multiparentSNP(keepHets = FALSE)
#' mapped <- mpcrossMapped(crossSNP, map = simulatedFourParentMap)
#' probabilities <- computeGenotypeProbabilities(mapped, error = 0.05)
#' probabilityData <- probabilityData(probabilities)
#' probabilityData[1:5, 1:5]
#' transposeProbabilities(probabilities)[1:5,1:5]
#' @export
transposeProbabilities <- function(inputObject)
{
	if(inherits(inputObject, "mpcrossMapped"))
	{
		if(length(inputObject@geneticData) == 1)
		{
			geneticData <- inputObject@geneticData[[1]]
		}
		else stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object")
	}
	else if(inherits(inputObject, "geneticData"))
	{
		if(is.null(inputObject@probabilities))
		{
			stop("Input object had no probabilities")
		}
		geneticData <- inputObject
	}
	else
	{
		stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object")
	}
	probabilities <- geneticData@probabilities
	nFounders <- nFounders(geneticData)
	nFinals <- nLines(geneticData)
	nProbabilitiesPositions <- length(unlist(probabilities@map))
	key <- probabilities@key
	nAlleles <- max(key[,3])
	if(geneticData@pedigree@selfing == "infinite" && (nrow(probabilities@data) != nFounders * nFinals || ncol(probabilities@data) != nProbabilitiesPositions))
	{
		stop("Slot probabilities@data had the wrong dimensions")
	}
	else if(geneticData@pedigree@selfing == "finite" && (nrow(probabilities@data) != nAlleles * nFinals || ncol(probabilities@data) != nProbabilitiesPositions))
	{
		stop("Slot probabilities@data had the wrong dimensions")
	}
	return(.Call("transposeProbabilities", geneticData, PACKAGE="mpMap2"))
}
