#' @export
as.mpInterval <- function(object, type = "mpMarker")
{
	if(!isS4(object) || !is(object, "mpcrossMapped"))
	{
		stop("Input object must be an S4 object of class mpcrossMapped")
	}
	if(type != "mpMarker")
	{
		stop("The only value for type is currently \"mpMarker\". Value \"mpInterval\" will be implemented later")
	}
	for(dataset in object@geneticData)
	{
		if(dataset@pedigree@selfing != "infinite")
		{
			stop("Only data-sets with infinite generations of selfing can be converted to objects of class mpInterval")
		}
		if(is.null(dataset@probabilities))
		{
			stop("At least one contained data-set did not have genotype probabilities. Please run computeGenotypeProbabilities first")
		}
	}
	nChromosomes <- length(object@map)
	chromosomeNames <- names(object@map)
	results <- list()
	flattenedMap <- unlist(lapply(object@map, names))
	for(datasetCounter in 1:length(object@geneticData))
	{
		geneticData <- object@geneticData[[datasetCounter]]
		currentDataFounders <- nFounders(geneticData)
		wgaimObject <- list()
		wgaimObject$pheno <- data.frame(id = rownames(finals(geneticData)))
		wgaimObject$geno <- vector("list", nChromosomes)
		names(wgaimObject$geno) <- chromosomeNames
		wgaimObject$nfounders <- nFounders(geneticData)
		wgaimObject$founders <- rownames(founders(geneticData))
		transformed <- .Call("transformForMPWGAIM", geneticData, PACKAGE="mpMap2")
		flattenedProbabilitiesMap <- unlist(lapply(geneticData@probabilities@map, names))
		for(chromosome in 1:length(object@map))
		{
			currentChrMap <- object@map[[chromosome]]
			currentChrProbabilitiesMap <- geneticData@probabilities@map[[chromosome]]
			wgaimObject$geno[[chromosome]]$map <- wgaimObject$geno[[chromosome]]$dist <- currentChrMap
			wgaimObject$geno[[chromosome]]$data <- transformed$finals[,names(currentChrMap)]
			wgaimObject$geno[[chromosome]]$founders <- transformed$founders[,names(currentChrMap)]
			markerIndices <- match(names(currentChrProbabilitiesMap), flattenedProbabilitiesMap)
			start <- min(markerIndices)
			end <- max(markerIndices)
			wgaimObject$geno[[chromosome]]$imputed.data <- transformed$probabilities[, (start*currentDataFounders - (currentDataFounders - 1)):(currentDataFounders*end)]
			class(wgaimObject$geno[[chromosome]]) <- "A"
		}
		class(wgaimObject) <- c("mpInterval", "cross", "interval")
		results[[length(results)+1]] <- wgaimObject
	}
	return(results)
}
