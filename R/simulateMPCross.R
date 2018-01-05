#' Simulate data from multi-parent designs
#' 
#' Data is simulated according to a pedigree, map and QTL model
#' @export
#' @param map Linkage map with which to generate data. See \code{\link[qtl]{sim.map}}
#' @param pedigree Pedigree for a multi-parent cross. 
#' @param mapFunction Map function used to convert distances to recombination fractions
#' @param seed Random seed to use. 
#' @return Object of class \code{mpcross}. 
#' @export
simulateMPCross <- function(map, pedigree, mapFunction, seed)
{
	isDetailedPedigreeArgument(pedigree)
	isMapArgument(map)
	if(!missing(seed))
	{
		nonNegativeIntegerArgument(seed)
		set.seed(seed)
	}

	#Treat markers as one big chromosome, with some recombination fractions of 0.5 to mark the gaps between actual chromosomes
	adjacentRecombination <- vector(mode="numeric")
	for (i in 1:length(map))
	{
		if(length(diff(map[[i]])) > 0)
		{
			map[[i]] <- sort(map[[i]])
			adjacentRecombination <- c(adjacentRecombination, sapply(diff(map[[i]]), mapFunction), 0.5)
		}
		else
		{
			adjacentRecombination <- c(adjacentRecombination, 0.5)
		}
	}
	nMarkers <- length(adjacentRecombination)
	markerNames <- unlist(lapply(map, names))
	#Remove last value of 0.5
	adjacentRecombination <- adjacentRecombination[-length(adjacentRecombination)]
	genotypes <- .Call("generateGenotypes", adjacentRecombination, markerNames, pedigree, PACKAGE="mpMap2")
	
	hetData <- fullHetData(map, nFounders(pedigree))

	#At this point we're still tracking both alleles. For the founders the alleles are the same, so we drop the second half
	founders <- genotypes[1:nFounders(pedigree),1:nMarkers,drop=FALSE]
	#For the finals we have to combine the two
	finals <- genotypes[which(pedigree@observed),,drop=FALSE]
	finalsRowNames <- rownames(finals)
	finalsColNames <- colnames(finals)[1:nMarkers]
	#Combine the pairs of alleles according to the hetData object
	finals <- combineGenotypes(finals, hetData)
	colnames(finals) <- finalsColNames
	rownames(finals) <- finalsRowNames
	return(new("mpcross", geneticData = new("geneticDataList", list(new("geneticData", founders = founders, finals = finals, hetData = hetData, pedigree = pedigree)))))
}
