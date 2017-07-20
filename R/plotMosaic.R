#' @export
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
		imputed <- imputed[lines, ]
	}
	dataMatrix <- imputed@data
	founders <- nFounders(geneticData)
	dataMatrix[dataMatrix > founders+1] <- founders+1
	heatmap_2(dataMatrix, scale = "none", col=brewer.pal(founders+1, "Spectral"), Colv = NA, Rowv = NA, do.dendro=c(FALSE,FALSE), ...)
}
