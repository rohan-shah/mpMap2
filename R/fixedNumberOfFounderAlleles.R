setClass("fixedNumberOfFounderAlleles", slots = list(alleles = "integer"))
#' @include geneticData-class.R
#' @include mpcross-class.R
#' @title Convert fully informative experiment to one with a fixed number of alleles per marker
#' @description Convert a fully informative experiment to one with a fixed number of alleles per marker
#' @details By default, simulated data is fully informative, so every founder carries its own allele, and all heterozygotes are distinguishable. 
#' 
#' This function takes in a fully informative experiment, and changes every marker so that it has a fixed number of founder alleles. Heterozygotes are also changed, so every combination of different alleles is still distinguishable. 
#' @param alleles Number of alleles for each marker
#' @return An object of internal class \code{fixedNumberOfFounderAlleles} suitable for application to an object of class \code{mpcross} using the addition operation. 
#' @examples 
#' data(simulatedFourParentData)
#' founders(simulatedFourParentData)[, 1:10]
#' altered <- simulatedFourParentData + fixedNumberOfFounderAlleles(3)
#' founders(altered)[, 1:10]
#' @export
fixedNumberOfFounderAlleles <- function(alleles)
{
	if(!is.numeric(alleles) || length(alleles) != 1) stop("Input alleles must be a single integer")
	return(new("fixedNumberOfFounderAlleles", alleles = as.integer(alleles)))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("geneticData", "fixedNumberOfFounderAlleles"), definition = function(e1, e2)
{
	nFounders <- nFounders(e1)
	expectedAlleleCount <- nFounders*(nFounders+1)/2
	if(e2@alleles >= nFounders)
	{
		stop("New number of founder alleles must be strictly less than the number of founders")
	}
	correctAlleles <- unlist(lapply(e1@hetData, function(x)
		{
			unique <- unique(x[,3])
			length(unique) == expectedAlleleCount && all(sort(unique) == 1:expectedAlleleCount)
		}))
	if(any(!correctAlleles)) stop("Can only apply fixedNumberOfFounderAlleles to fully informative designs where the encodings are consecutive numbers starting at 1")
	markers <- nMarkers(e1)
	copies <- as.integer(floor(nFounders / e2@alleles))
	leftOver <- nFounders - copies*e2@alleles
	copiesTable <- rep(1:e2@alleles, times = copies)

	totalNewHetDataRows <- e2@alleles * e2@alleles
	totalNewHets <- (totalNewHetDataRows - e2@alleles)/2

	newHetData <- matrix(0L, nrow = totalNewHetDataRows, ncol = 3)
	newHetData[1:e2@alleles,] <- rep(1:e2@alleles, times = 3)
	newHetData[(e2@alleles+1):(e2@alleles + totalNewHets),3] <- (e2@alleles+1):(e2@alleles + totalNewHets)
	newHetData[(e2@alleles+1):(e2@alleles + totalNewHets), 1:2] <- t(combn(1:e2@alleles, 2))
	newHetData[(e2@alleles+1+totalNewHets):(e2@alleles + 2*totalNewHets), 1] <- newHetData[(e2@alleles+1):(e2@alleles + totalNewHets), 2]
	newHetData[(e2@alleles+1+totalNewHets):(e2@alleles + 2*totalNewHets), 2] <- newHetData[(e2@alleles+1):(e2@alleles + totalNewHets), 1]
	newHetData[(e2@alleles+1+totalNewHets):(e2@alleles + 2*totalNewHets), 3] <- (e2@alleles+1):(e2@alleles + totalNewHets)
	sapply(1:markers, function(x)
		{
			table <- c(copiesTable, sample(1:e2@alleles, leftOver, replace=FALSE))
			table <- sample(table)
			newFounders <- table[e1@founders[,x]]
			e1@founders[,x] <<- newFounders
			oldHetData <- e1@hetData[[x]]
			e1@hetData[[x]] <<- newHetData

			finalsTable <- sapply(1:expectedAlleleCount, function(y)
				{
					relevantOldRow <- which(oldHetData[,3] == y)[1]
					newFounderEncodings <- c(table[oldHetData[relevantOldRow,1]], table[oldHetData[relevantOldRow,2]])
					relevantNewRow <- which(newHetData[,1] == newFounderEncodings[1] & newHetData[,2] == newFounderEncodings[2])[1]
					return(newHetData[relevantNewRow,3])
				})
			e1@finals[,x] <<- finalsTable[e1@finals[,x]]
		})
	return(e1)
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "fixedNumberOfFounderAlleles"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Changing marker patterns will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	if(length(e1@geneticData) > 1)
	{
		stop("Attempting to change an object containing multiple data sets. Please change each dataset individually")
	}
	e1@geneticData[[1]] <- e1@geneticData[[1]]+e2
	return(e1)
})
