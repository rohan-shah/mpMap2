#' @export
transposeProbabilities <- function(probabilities)
{
	if(inherits(probabilities, "mpcrossMapped"))
	{
		mpcross <- probabilities
		if(length(mpcross@geneticData) == 1)
		{
			probabilities <- mpcross@geneticData[[1]]
		}
		else stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object, or a probabilities object")
	}
	if(inherits(probabilities, "geneticData"))
	{
		if(!is.null(probabilities@probabilities))
		{
			stop("Input object had no probabilities")
		}
		probabilities <- probabilities@probabilities
	}
	if(!inherits(probabilities, "probabilities"))
	{
		stop("Please input an mpcrossMapped object with a single experiment, or a geneticData object, or a probabilities object")
	}
	nFounders <- length(unique(probabilities@key[,1]))
	nFinals <- nrow(probabilities@data) / nAlleles
	nProbabilitiesPositions <- length(unlist(probabilities@map))
	key <- probabilities@key
	nAlleles <- max(key[,3])
	if(nrow(probabilities@data) != nAlleles * nFinals || ncol(probabilities@data) != nProbabilitiesPositions)
	{
		stop("probabilities@data had the wrong dimensions")
	}
	return(.Call("transposeProbabilities", probabilities, PACKAGE="mpMap2"))
}
