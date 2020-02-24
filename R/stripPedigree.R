#' @title Strip pedigree of unneccessary lines
#' @description Strip pedigree of lines that make no genetic contribution to the specified set of lines.
#' @details Pedigrees for structured experiments can be messy. Often they include lines that make no genetic contribution to the lines that were finally genotyped. When it comes to visualising the structure of the pedigree of the final population, these unneccessary extra lines can make it difficult to see the structure. This function takes in a pedigree and a list of genetic lines, and returns a subpedigree that contains only those lines that make a genetic contribution to the named lines.
#'
#' This function relies on the use of the Boost C++ libraries, and may not be available in every distributed version of mpMap2. If this function is unavailable, the function will return \code{NULL}. 
#' @param pedigree The initial pedigree, which may contain some unneccessary extra genetic lines
#' @param finalLines The list of lines of interest. Lines in the pedigree which do not make a genetic contribution to the lines in \code{finalLines} will be removed. 
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
stripPedigree <- function(pedigree, finalLines)
{
	isPedigreeArgument(pedigree)
 	if(!is.character(finalLines) || length(finalLines) == 0)
	{
		stop("Input finalLines must be a character vector")
	}
	success <- FALSE
	try(
	{
		newPedigree <- .Call("stripPedigree", pedigree, finalLines, PACKAGE="mpMap2")
		success <- TRUE
	}, silent=TRUE)
	if(success) return(newPedigree)
	return(NULL)
}
