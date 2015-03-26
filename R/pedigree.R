lineByName <- function(pedigree, name)
{
	isPedigreeArgument(pedigree)
	lineIndex <- match(name, pedigree@lineNames)
	if(is.na(lineIndex))
	{
		stop(paste0("Could not find line named ", name))
	}
	return(c(lineName = name, mother = pedigree@lineNames[pedigree@mother[lineIndex]], father = pedigree@lineNames[pedigree@father[lineIndex]]))
}
setMethod(f = "length", signature = "pedigree", definition = function(x)
{
	return(length(x@lineNames))
})
setMethod(f = "print", signature = "pedigree", definition = function(x)
{
	cat("This is a pedigree object containing ", length(x), " lines, of which ", nFounders(x), " are founders and ", sum(x@observed), " are observed\n")
})