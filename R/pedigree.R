linesByNames <- function(pedigree, names)
{
	isPedigreeArgument(pedigree)
	lineIndices <- match(names, pedigree@lineNames)
	if(any(is.na(lineIndices)))
	{
		if(sum(is.na(lineIndices)) == 1) stop(paste0("Could not find line named ", names[which(is.na(lineIndices))]))
		stop(paste0("Could not find lines named ", do.call(paste, as.list(c(names[which(is.na(lineIndices))], sep = ", ")))))
	}
	return(cbind(lineName = names, mother = pedigree@lineNames[pedigree@mother[lineIndices]], father = pedigree@lineNames[pedigree@father[lineIndices]]))
}
setMethod(f = "print", signature = "pedigree", definition = function(x)
{
	cat("This is a pedigree object containing ", length(x), " lines\n")
})
setMethod(f = "print", signature = "detailedPedigree", definition = function(x)
{
	cat("This is a pedigree object containing ", length(x@lineNames), " lines, of which ", nFounders(x), " are founders and ", nFinals(x), " are observed\n")
})