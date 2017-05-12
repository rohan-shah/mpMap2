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
