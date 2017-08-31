#' @export
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
#' @include pedigree-class.R
#' @include detailedPedigree-class.R
setMethod(f = "print", signature = "pedigree", definition = function(x)
{
	cat("This is a pedigree object containing ", length(x), " lines\n")
})
setMethod(f = "print", signature = "detailedPedigree", definition = function(x)
{
	cat("This is a pedigree object containing ", length(x@lineNames), " lines, of which ", nFounders(x), " are founders and ", nLines(x), " are observed\n")
})
#' @export
pedigree <- function(lineNames, mother, father, selfing, warnImproperFunnels = TRUE)
{
	if(!is.numeric(mother))
	{
		stop("Input mother must be a numeric vector")
	}
	if(!is.numeric(father))
	{
		stop("Input father must be a numeric vector")
	}
	if(!is.character(lineNames))
	{
		stop("Input lineNames must be a character vector")
	}
	if(length(mother) != length(father) || length(father) != length(lineNames))
	{
		stop("Inputs mother, father and lineNames must have the same length")
	}
	if(selfing != "infinite" && selfing != "finite")
	{
		stop("Input selfing must be either \"finite\" or \"infinite\"")
	}
	mode(mother) <- "integer"
	mode(father) <- "integer"
	pedigree <- new("pedigree", mother = mother, father = father, lineNames = lineNames, selfing = selfing, warnImproperFunnels = warnImproperFunnels)
	return(pedigree)
}
