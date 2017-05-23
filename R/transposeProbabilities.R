#' @export
transposeProbabilities <- function(geneticData)
{
	if(inherits(geneticData, "mpcrossMapped"))
	{
		mpcross <- geneticData
		if(length(mpcross@geneticData) == 1)
		{
			geneticData <- mpcross@geneticData[[1]]
		}
		else stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object")
	}
	if(!inherits(geneticData, "geneticData"))
	{
		stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object")
	}
	nFounders <- nrow(founders(geneticData))
	nFinals <- nrow(finals(geneticData))
	nProbabilitiesPositions <- length(unlist(geneticData@probabilities@map))
	key <- geneticData@probabilities@key
	nAlleles <- max(key[,3])
	if(nrow(geneticData@probabilities@data) != nAlleles * nFinals || ncol(geneticData@probabilities@data) != nProbabilitiesPositions)
	{
		stop("geneticData@probabilities@data had the wrong dimensions")
	}
	return(.Call("transposeProbabilities", geneticData, PACKAGE="mpMap2"))
}
