#' @export
generateGridPositions <- function(spacing)
{
	retFunction <- function(object)
	{
		result <- lapply(as.list(names(object@map)), function(chrName)
			{
				x <- object@map[[chrName]]
				range <- range(x)
				positions <- seq(range[1], range[2], by = spacing)
				names(positions) <- paste0("Chr", chrName, "Loc", 1:length(positions))
				positions
			})
		names(result) <- names(object@map)
		result
	}
	return(retFunction)
}
#' @export
imputeFounders <- function(mpcrossMapped, homozygoteMissingProb = 1, heterozygoteMissingProb = 1, errorProb = 0, extraPositions = list(), showProgress = FALSE)
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
		results <- .Call("imputeFounders", mpcrossMapped@geneticData[[i]], mpcrossMapped@map, homozygoteMissingProb, heterozygoteMissingProb, errorProb, extraPositions, showProgress, PACKAGE="mpMap2")
		resultsMatrix <- results$data
		errors <- NULL
		if(errorProb != 0)
		{
			errors <- results$errors
			rownames(errors) <- rownames(resultsMatrix)
			colnames(errors) <- colnames(resultsMatrix)
		}
		class(results$map) <- "map"
		names(results$map) <- names(map)
		mpcrossMapped@geneticData[[i]]@imputed <- new("imputed", data = resultsMatrix, key = results$key, map = results$map, errors = errors)
	}
	return(mpcrossMapped)
}
