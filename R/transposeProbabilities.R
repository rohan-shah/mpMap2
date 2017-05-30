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
		stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object, or a probabilities object")
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
