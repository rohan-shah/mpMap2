#' @export
imputeFounders <- function(mpcrossMapped, homozygoteMissingProb = 1, heterozygoteMissingProb = 1, errorProb = 0, extraPositions = list())
{
	isNewMpcrossMappedArgument(mpcrossMapped)
	if(homozygoteMissingProb < 0 || homozygoteMissingProb > 1)
	{
		stop("Input homozygoteMissingProb must be a value between 0 and 1")
	}
	if(heterozygoteMissingProb < 0 || heterozygoteMissingProb > 1)
	{
		stop("Input heterozygoteMissingProb must be a value between 0 and 1")
	}
	if(errorProb < 0 || errorProb >= 1)
	{
		stop("Input errorProb must be non-negative and smaller than 1")
	}
	map <- mpcrossMapped@map
	if(!all(names(extraPositions) %in% names(map)))
	{
		stop("Input extraPositions must be a list, with entries named after chromosomes")
	}
	allMarkerNames <- unlist(lapply(map, names))
	for(chromosome in names(extraPositions))
	{
		extraChr <- extraPositions[[chromosome]]
		if(any(names(extraChr) %in% allMarkerNames))
		{
			stop("Extra locations in extraPositions cannot be named after markers")
		}
		if(!is.numeric(extraChr))
		{
			stop("Input extraPositions must be a list, with entries which are numeric vectors")
		}
		if(is.null(names(extraChr)))
		{
			stop("Vectors in input extraPositions must have names")
		}
		extraPositions[[chromosome]] <- sort(extraChr)
	}
	for(i in 1:length(mpcrossMapped@geneticData))
	{
		results <- .Call("imputeFounders", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, heterozygoteMissingProb, errorProb, extraPositions, PACKAGE="mpMap2")
		resultsMatrix <- results$data
		mpcrossMapped@geneticData[[i]]@imputed <- new("imputed", data = resultsMatrix, key = results$key)
	}
	return(mpcrossMapped)
}
