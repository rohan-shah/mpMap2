purdyToPedigree <- function(lineNames, purdy, selfing, warnImproperFunnels)
{
	if(!is.character(lineNames) || !is.character(purdy))
	{
		stop("Inputs lineNames and purdy must be character vectors")
	}
	if(length(lineNames) != length(purdy))
	{
		stop("Inputs lineNames and purdy must have the same lengths")
	}
	parsed <- .Call("parsePurdy", lineNames, purdy, PACKAGE="mpMap2")
	lineNames <- parsed[,1]
	mother <- sapply(parsed[,2], function(x)
	{
		if(x == "") return(0)
		return(match(x, lineNames))
	})
	father <- sapply(parsed[,3], function(x)
	{
		if(x == "") return(0)
		return(match(x, lineNames))
	})
	names(father) <- names(mother) <- NULL
	pedigree <- pedigree(lineNames = lineNames, mother = mother, father = father, selfing = selfing, warnImproperFunnels = warnImproperFunnels)
	return(pedigree)
}
