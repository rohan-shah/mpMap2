reorderPedigree <- function(mother, father, lineNames, warnImproperFunnels, selfing)
{
	success <- FALSE
	try(
	{
		result <- .Call("reorderPedigree", mother, father, lineNames, PACKAGE="mpMap2")
		success <- TRUE
	}, silent=TRUE)
	#There will be an error if package was compiled without boost support
	if(success)
	{
		return(new("pedigree", lineNames = result$lineNames, mother = result$mother, father = result$father, selfing = selfing, warnImproperFunnels = warnImproperFunnels))
	}
	else
	{
		return(NULL)
	}
}
