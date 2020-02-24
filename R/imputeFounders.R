#' @export
#' @title Impute underlying genotypes
#'
#' @description
#' Impute the most likely sequence of underlying genotypes, using the Viterbi algorithm
#'
#' @param mpcrossMapped An object containing genetic data and a genetic map
#' @param homozygoteMissingProb The probability with which homozygous genotypes are observed as missing.
#' @param heterozygoteMissingProb The probability with which heterozygous genotypes are observed as missing.
#' @param errorProb The probability of a genotyping error. 
#' @param extraPositions Extra genetic positions at which to perform imputation. 
#' @param showProgress If this paramater is \code{TRUE}, a progress bar is produced. 
#' @return An object of class \code{mpcrossMapped}, containing all the information in the input object, and also including imputed IBD genotypes. 
#' This function uses the Viterbi algorithm to calculate the most likely sequence of underlying genotypes, given observed genetic data. The parameters for the algorithm are a homozygous mising rate, a heterozygous missing rate, and an error probability. 
#' 
#' The two missing rates are intended to allow long strings of missing values to be imputed as heterozygotes, in the case that heterozygous genotypes are observed as missing much more often than homozygotes. Only the ratio of these two parameters is relevant, which is why the default values of 1 are acceptable. These default values really mean that the missing rates are equal. 
#' 
#' The parameter \code{extraPositions} specifies the genetic positions at which imputation should be performed. This can be either a list, or a function such as \code{generateGridPositions} \code{generateIntervalMidPoints}. If a function is input, this function is applied to the input genetic map, to determine the extra genetic locations. If a list is input, the names of the list entries should be chromosome names, and the entry for each chromosome should be a named vector. We give an example of the list format in the examples section at the bottom of this page. 
#' 
#' One subtlety when using extra genetic positions is that specifying such positions can change the results of the imputation process. This is undesirable, but does not represent a bug in the implementation. The Hidden Markov Model (HMM) used to model the genotypes is not exact, although it is a highly accurate approximation. As it is an approximation, it fails to satisfy the condition 
#' \deqn{P^{s+t} = P^t P^s}
#' This property (a stochastic semigroup property) fails to hold because the HMM is only an approximation. As a result, adding extra genetic positions can change the results of the imputation. We emphasise that this is possible only when there are number of underlying sequences which are almost equally likely, and even then this problem occurs rarely. However, this problem becomes obvious when large simulation studies are performed. 
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
