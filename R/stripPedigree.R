#' @export
stripPedigree <- function(pedigree, finalLines)
{
	isPedigreeArgument(pedigree)
 	if(!is.character(finalLines) || length(finalLines) == 0)
	{
		stop("Input finalLines must be a character vector")
	}
	newPedigree <- .Call("stripPedigree", pedigree, finalLines, PACKAGE="mpMap2")
	return(newPedigree)
}
