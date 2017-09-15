computeGenotypeProbabilitiesInternal <- function(geneticData, map, homozygoteMissingProb, heterozygoteMissingProb, errorProb, extraPositions)
{
	results <- .Call("computeGenotypeProbabilities", geneticData, map, homozygoteMissingProb, heterozygoteMissingProb, errorProb, extraPositions, PACKAGE="mpMap2")
	resultsMatrix <- results$data
	founderNames <- rownames(geneticData@founders)
	if(geneticData@pedigree@selfing == "infinite")
	{
		rownames(resultsMatrix) <- unlist(lapply(rownames(geneticData@finals), function(lineName) paste0(lineName, " - ", founderNames)))
	}
	else
	{
		nAlleles <- nrow(resultsMatrix) / nrow(geneticData@finals)
		rownames(resultsMatrix) <- unlist(lapply(rownames(geneticData@finals), function(lineName) paste0(lineName, " - ", 1:nAlleles)))
	}
	class(results$map) <- "map"
	names(results$map) <- names(map)
	return(new("probabilities", data = resultsMatrix, key = results$key, map = results$map))
}
#' @export
computeGenotypeProbabilities <- function(mpcrossMapped, homozygoteMissingProb = 1, heterozygoteMissingProb = 1, errorProb = 0, extraPositions = list())
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
	map <- mpcrossMapped@map
	allMarkerNames <- unlist(lapply(map, names))

	#Input extraPositions can be a list or a function
	if(class(extraPositions) != "list" && class(extraPositions) != "function")
	{
		stop("Input extraPositions must be a list or a function that generates a list")
	}
	if(class(extraPositions) == "function")
	{
		extraPositions <- extraPositions(mpcrossMapped)
	}
	if(!all(names(extraPositions) %in% names(map)))
	{
		stop("Input extraPositions must be a list, with entries named after chromosomes")
	}
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
		mpcrossMapped@geneticData[[i]]@probabilities <- computeGenotypeProbabilitiesInternal(mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, heterozygoteMissingProb, errorProb, extraPositions)
	}
	return(mpcrossMapped)
}
