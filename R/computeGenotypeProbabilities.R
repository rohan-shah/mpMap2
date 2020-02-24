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
#' @title Compute IBD genotype probabilities
#' @description Compute IBD genotype probabilities
#' @details This function computes the IBD genotype probabilities using a Hidden Marker Model (HMM) and the forward-backward algorithm. The HMM model is only an approximation to the underlying genetics, but it is a very good one. 
#' 
#' There are a number of parameters to this model. \code{homozygoteMissingProb} gives the "probability" that a marker homozygote will be marked a missing. \code{heterozygoteMissingProb} gives the "probability" that a marker heterozygote will be marked as missing. We say "probability" because really the important thing is the difference these two parameters, not the values themselves. If they are equal then a missing marker genotype contains no information. If code{heterozygoteMissingProb} is relatively larger than \code{homozygoteMissingProb}, then missing marker genotypes suggests that the underlying genotype is heterozygous, provided enough missing marker values occur sequentially. 
#' 
#' The key reason for introducing these paramters is that if \code{heterozygoteMissingProb} is relatively larger, then a dataset with no observed marker heterozygotes can still be used to estimate positions of underlying heterozygous genotypes, provided that heterozygous genotypes lead to consecutive missing marker genotype values.
#' 
#' The \code{errorProb} parameter gives the probability that a marker genotype is actually incorrect. In this case, it is assumed that the correct value for this marker genotype is random and uniformly distributed. This is different from assuming that the underlying genotype itself is random. If \code{errorProb} is zero, then it is not possible to have co-located markers with inconsistent genotypes, and if this occurs an error will be generated. \code{jitterMap} can be used to avoid this, but setting \code{errorProb} to some non-zero value is a much better solution. 
#' 
#' It is also possible to generate IBD probabilities at non-marker positions. These extra positions are specified by the \code{extraPositions} options, which can be specified two ways. The first is by specifying a list with name entries, where the names correspond to chromosomes. Each named entry should be a named vector, with names corresponding to the names of the positions, and values corresponding to the positions in cM on that chromosome. 
#'
#' The second possibility is to specify a function, which will be applied to the input object of class \code{mpcrossMapped} to generate the extra positions. Two helper options are provided for this - \code{\link{generateGridPositions}} and \code{link{generateIntervalMidPoints}}. 
#' 
#' @param mpcrossMapped An object of class \code{mpcrossMapped}, for which to compute the IBD genotype probabilities.
#' @param homozygoteMissingProb The "probability" that a marker genotype that is truly homozygous will be marked as missing.
#' @param heterozygoteMissingProb The "probability" that a marker genotype that is truly heterozygous will be marked as missing.
#' @param errorProb The probability that a marker genotype is incorrect.
#' @param extraPositions The extra positions at which to compute the IBD genotype probabilities. May be either a list with named components corresponding to chromosomes (simialr to a map) or a function which will be applied to the input object to generate the extra positions.
#' @return An object of class \code{mpcrossMapped} containing all information in the input object, and also estimated IDB probabilities.
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
