#' @title Extract pedigree by names
#' @description Extract part of pedigree in human-readable format
#' @details Pedigrees in mpMap2 are stored using indices for maternal and paternal lines, which is not a human-readable format. This function takes in a pedigree, and returns a human-readable subset.
#' @param pedigree An object of class \code{pedigree}
#' @param names The names of the lines for which to extract the pedigree
#' @return A matrix giving the genetic lines and their parents, by line name. 
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
#' @title Create a pedigree object
#' @description Create a pedigree object
#' @details This function creates a pedigree object from parts. All lines are assumed to have an index, starting at 1 for the first line. Values at index of the various inputs 1 all relate to the first line, values at index 2 all relate to the second line, etc. 
#' 
#' Input \code{lineNames} assigns a name to every line. Input \code{mother} gives the index of a mother line, where a value of 0 indicates that a line is a founder of the population (and therefore inbred). Input \code{father} gives the index of a father line, where a value of 0 indicates that a line is a founder of the population (and therefore inbred). Input \code{selfing} must be \code{"finite"} or \code{"infinite"}. A value of infinite means that the number of generations of selfing for this pedigree will be assumed to be infinite. A value of \code{"finite"} means that the number of generations of selfing will be computed from the pedigree, for every line.  
#' @param lineNames The names of the genetic lines
#' @param mother The index of the maternal line
#' @param father The index of the paternal line
#' @param selfing Should the number of generations of selfing be taken from the pedigree (\code{"finite"}), or should selfing be assumed to be infinite (\code{"infinite"})? 
#' @param warnImproperFunnels Should a warning be generated in subsequent computations using this pedigree, if there are lines which do not contain all founding lines as ancestors?
#' @return An object of class \code{pedigree} representing the inputs.
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
