#' @export
#' @title Test for distortion using IBD genotype probabilities
#' @description Test for distortion using IBD genotype probabilities
#' @details In real experiments, genetic inheritance may not follow the expected model. This function tests for deviations from expected inheritance by using the genetic composition of the population at individual positions, as measured by the IBD genotype probabilities. 
#' 
#' At a particular point, the mean for each founder allele of the IBD genotype probabilities for each founder allele are summed across the population. The average is taken, and this is then compared with the proportion expected to be inherited from that founder, under standard models of genetic inheritance. 
#' 
#' The result is a matrix containing p-values, test-statistic values, and the L1 and L2 distances between the observed genetic proportions, and the expected genetic proportions. 
#' @param object An object of class \code{mpcrossMapped} which contains imputed IBD genotype data
#' @return A data.frame containing p-values and test-statistic values for each position at which there is IBD genotype probability data. 
testDistortion <- function(object)
{
	if(inherits(object, "mpcross") && length(object@geneticData) == 1)
	{
		object <- object@geneticData[[1]]
	}
	if(!isS4(object) || !inherits(object, "geneticData"))
	{
		stop("Input object must be an S4 object of class geneticData")
	}
	result <- .Call("testDistortion", object, PACKAGE="mpMap2")
	pValues <- 1 - pchisq(result$testStatistics, result$classes - 1)
	result <- data.frame(pValue = pValues, statistic = result$testStatistics, L1 = result$L1, L2 = result$L2)
	rownames(result) <- markers(object)
	return(result)
}
