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