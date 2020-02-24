#' @export
#' @title Plot estimated genetic composition of lines
#' @description Plot estimated genetic composition of lines
#' @details This function produces a heatmap showing the genetic composition of lines, as measured by the imputed IBD genotypes. Rows correspond to genetic lines, columns correspond to genetic positions, and colours indicate founder alleles. All heterozygotes are marked in the same colour, otherwise there are generally too many colours to be useful. 
#' @param inputObject An object of class \code{mpcrossMapped} containing imputed IBD genotypes 
#' @param chromosomes Chromosomes to plot
#' @param positions Genetic positions to plot
#' @param lines Genetic lines to plot
#' @param ... Extra inputs to \code{heatmap_2}
#' @return None
plotMosaic <- function(inputObject, chromosomes, positions, lines, ...)
{
	if(inherits(inputObject, "mpcrossMapped"))
	{
		if(length(inputObject@geneticData) != 1)
		{
			stop("Input object contained more than one data set")
		}
		geneticData <- inputObject@geneticData[[1]]
	}
	else if(inherits(inputObject, "geneticData"))
	{
		geneticData <- inputObject
	}
	else
	{
		stop("Input object must have class mpcrossMapped or geneticData")
	}
	if(!missing(positions) && !missing(chromosomes))
	{
		stop("Cannot use inputs positions and chromosomes together")
	}
	if(!missing(chromosomes))
	{
		subsettedMap <- inputObject@map[names(inputObject@map) %in% chromosomes]
		if(length(subsettedMap) == 0 || any(!(chromosomes %in% names(inputObject@map)))) stop("Please enter valid chromosome names")
		positions <- unlist(lapply(subsettedMap, names))
	}
	if(!missing(positions))
	{
		imputed <- subset(geneticData@imputed, positions = positions)
	}
	else
	{
		imputed <- geneticData@imputed
	}
	if(!missing(lines))
	{
		imputed <- subset(imputed, lines = lines)
	}
	dataMatrix <- imputed@data
	founders <- nFounders(geneticData)
	dataMatrix[dataMatrix > founders+1] <- founders+1
	if(nrow(dataMatrix) == 1)
	{
		heatmap_2(rbind(dataMatrix, dataMatrix), scale = "none", col=brewer.pal(founders+1, "Spectral"), Colv = NA, Rowv = NA, do.dendro=c(FALSE,FALSE), ...)
	}
	else
	{
		heatmap_2(dataMatrix, scale = "none", col=brewer.pal(founders+1, "Spectral"), Colv = NA, Rowv = NA, do.dendro=c(FALSE,FALSE), ...)
	}
}
