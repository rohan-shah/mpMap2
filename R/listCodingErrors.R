#' @title Generate a list of encoding errors
#' @description Generate a list of encoding errors from genetic data
#' @export
#' @details Given genetic data matrices for the founding lines and the final lines of a population, and information about the encoding of marker heterozygotes, generate a list of errors. 
#' These errors include observed values which don't correspond to a known combination of marker alleles, missing values in the genetic data for the founding lines, etc. 
#' 
#' The results of this function allow human-readable lists of errors to be generated, or errors to be automatically fixed (if the errors are sufficiently simple). 
#' @param founders Genetic data for the founding lines of the population
#' @param finals Genetic data for the final lines of the population
#' @param hetData Data about the encoding of marker heterozygotes
#' @return List with the following entries:
#' 	\describe{
#' 		\item{finals}{Markers with an invalid observed value.}
#'		\item{null}{Markers with a missing value for a founding line, for which the are observations for at least one genetic line.}
#' 		\item{missingHetData}{Markers for which a homozygote did not have an encoding.}
#'		\item{invalidHetData}{Markers for which the heterozygote encoding data was invalid.}
#'	}
listCodingErrors <- function(founders, finals, hetData)
{
	errors <- .Call("listCodingErrors", founders, finals, hetData, PACKAGE="mpMap2")
	errors$finals <- errors$finals + 1
	errors$invalidHetData <- errors$invalidHetData + 1
	errors$missingHetData[,1] <- errors$missingHetData[,1] + 1
	errors$null <- errors$null + 1
	return(errors)
}
#' @title Generate a list of encoding errors assuming infinite selfing
#' @description Generate a list of encoding errors assuming infinite selfing 
#' @details Generate a list of encoding errors assuming infinite selfing. Given the infinite selfing assumption, no information about heterozygote encoding is required. 
#' @param founders Genetic data for the founding lines of the population
#' @param finals Genetic data for the final lines of the population
#' @return List with the following entries:
#' 	\describe{
#' 		\item{finals}{Markers with an invalid observed value.}
#'		\item{null}{Markers with a missing value for a founding line, for which the are observations for at least one genetic line.}
#' 		\item{missingHetData}{Markers for which a homozygote did not have an encoding.}
#'		\item{invalidHetData}{Markers for which the heterozygote encoding data was invalid.}
#'	}
#' @export
listCodingErrorsInfiniteSelfing <- function(founders, finals)
{
	newHetDataList <- lapply(as.list(1:ncol(founders)), function(x)
	{
		if(any(is.na(founders[,x])))
		{
			finals[,x] <<- NA
			return(matrix(0L, 0, 3))
		}
		else
		{
			uniqueAlleles <- unique(founders[, x])
			retVal <- cbind(uniqueAlleles, uniqueAlleles, uniqueAlleles)
			colnames(retVal) <- NULL
			return(retVal)
		}
	})
	names(newHetDataList) <- colnames(founders)

	newHetData <- new("hetData", newHetDataList)
	errors <- .Call("listCodingErrors", founders, finals, newHetData, PACKAGE="mpMap2")
	errors$finals <- errors$finals + 1
	errors$invalidHetData <- errors$invalidHetData + 1
	errors$missingHetData[,1] <- errors$missingHetData[,1] + 1
	errors$null <- errors$null + 1
	return(errors)
}
