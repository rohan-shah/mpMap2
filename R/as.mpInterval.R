#' @export
#' @title Convert mpcross object to MPWGAIM format
#' @description Convert an object of class \code{mpcrossMapped} to the format used by MPWGAIM. 
#' @details MPWGAIM is a package for performing QTL analysis using multi-parent populations. This function outputs a data object suitable for input to MPWGAIM. The output object can be in MPWGAIMs \code{mpMarker} or \code{mpInterval} formats. See the documentation of MPWGAIM for further information. 
#' @param object The object of class \code{mpcrossMapped} to convert
#' @param type The type of MPWGAIM object to output. Must be \code{"mpMarker"} or \code{"mpInterval"}
#' @param positions In the case of \code{mpMarker} format, the positions at which the IBD probabilities should be output. Must be either \code{"all"} (all positions for which IBD probabilities are available) or \code{"marker"} (only marker positions). 
#' @param homozygoteMissingProb Used as an input to \code{computeGenotypeProbabilitiesInternal}, if the IBD probabilities need to be calculated.
#' @param heterozygoteMissingProb Used as an input to \code{computeGenotypeProbabilitiesInternal}, if the IBD probabilities need to be calculated.
#' @param errorProb Used as an input to \code{computeGenotypeProbabilitiesInternal}, if the IBD probabilities need to be calculated.
#' @return An object of class \code{mpMarker} or \code{mpInterval}, which are formats specified by package mpwgaim. 
as.mpInterval <- function(object, type = "mpMarker", positions, homozygoteMissingProb, heterozygoteMissingProb, errorProb)
{
	if(!isS4(object) || !is(object, "mpcrossMapped"))
	{
		stop("Input object must be an S4 object of class mpcrossMapped")
	}
	if(type == "mpMarker")
	{
		if(!missing(homozygoteMissingProb))
		{
			warning("Input homozygoteMissingProb is ignored if type is \"mpInterval\"")
		}
		if(!missing(heterozygoteMissingProb))
		{
			warning("Input heterozygoteMissingProb is ignored if type is \"mpInterval\"")
		}
		if(!missing(errorProb))
		{
			warning("Input errorProb is ignored if type is \"mpInterval\"")
		}
		as.mpInterval.mpMarker(object = object, positions = positions)
	}
	else if(type == "mpInterval")
	{
		if(!missing(positions))
		{
			warning("Input positions is ignored if type is \"mpInterval\"")
		}
		as.mpInterval.mpInterval(object = object, homozygoteMissingProb = homozygoteMissingProb, heterozygoteMissingProb = heterozygoteMissingProb, errorProb = errorProb)
	}
	else stop("Input type must be either \"mpMarker\" or \"mpInterval\"")
}
as.mpInterval.mpInterval <- function(object, homozygoteMissingProb, heterozygoteMissingProb, errorProb)
{
	warnedProbabilities <- FALSE
	midPoints <- generateIntervalMidPoints(object)
	midPointNames <- unlist(lapply(midPoints, names))

	nChromosomes <- length(object@map)
	chromosomeNames <- names(object@map)
	results <- list()
	flattenedMap <- unlist(lapply(object@map, names))
	for(i in 1:length(object@geneticData))
	{
		if(object@geneticData[[i]]@pedigree@selfing != "infinite")
		{
			stop("Only data-sets with infinite generations of selfing can be converted to objects of class mpInterval")
		}
		if(is.null(object@geneticData[[i]]@probabilities))
		{
			if(!warnedProbabilities)
			{
				warning("At least one contained data-set did not have genotype probabilities. Probabilities will be computed. ")
				warnedProbabilities <- TRUE
			}
			if(missing(homozygoteMissingProb) || missing(heterozygoteMissingProb) || missing(errorProb))
			{
				stop("Inputs homozygoteMissingProb, heterozygoteMissingProb and errorProb must be specified, if genotype probabilities are to be computed")
			}
			object@geneticData[[i]]@probabilities <- computeGenotypeProbabilitiesInternal(geneticData = object@geneticData[[i]], map = object@map, homozygoteMissingProb = homozygoteMissingProb, heterozygoteMissingProb = heterozygoteMissingProb, errorProb = errorProb, extraPositions = midPoints)
		}
		else
		{
			if(!all(midPointNames %in% colnames(object@geneticData[[i]]@probabilities@data)))
			{
				stop("Unable to find probability data for the mid-points of every interval")
			}
			probabilitiesSubsettedMap <- lapply(object@geneticData[[i]]@probabilities@map, function(x) x[names(x) %in% midPointNames])
			if(!all.equal(probabilitiesSubsettedMap, midPoints, check.attributes = FALSE))
			{
				stop("Unable to find probability data for the mid-points of every interval - Some interval midpoints were not in the expected positions")
			}
		}
		geneticData <- object@geneticData[[i]]
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

			positionIndices <- match(names(midPoints[[chromosome]]), flattenedProbabilitiesMap)
			probabilityPositionIndices <- unlist(sapply(positionIndices, function(x) (x*currentDataFounders - (currentDataFounders - 1)):(x*currentDataFounders)))
			wgaimObject$geno[[chromosome]]$intval <- transformed$probabilities[,probabilityPositionIndices]

			currentChrUnique <- currentChrMap[!duplicated(currentChrMap)]
			#There is a class attribute that needs to be kept
			class(currentChrUnique) <- class(currentChrMap)
			wgaimObject$geno[[chromosome]]$map <- wgaimObject$geno[[chromosome]]$dist <- currentChrUnique
			wgaimObject$geno[[chromosome]]$data <- transformed$finals[,names(currentChrUnique)]
			wgaimObject$geno[[chromosome]]$founders <- transformed$founders[,names(currentChrUnique)]
		}
		class(wgaimObject) <- c("mpInterval", "cross", "interval")
		wgaimObject$gen.type <- "mpInterval"
		results[[length(results)+1]] <- wgaimObject
	}
	if(length(results) == 1) return(results[[1]])
	return(results)
}
as.mpInterval.mpMarker <- function(object, positions)
{
	if(missing(positions))
	{
		stop("Input positions must be either \"all\" or \"marker\"")
	}
	if(!(positions %in% c("all", "marker")))
	{
		stop("Input positions must be either \"all\" or \"marker\"")
	}
	if(positions == "all")
	{
		warning("Using positions = \"all\" will result in genetic data being replaced by NA's")
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

			#Determine whether we keep the genetic data. If we're using *all* the locations at which the probabilities were calculated, then there isn't necessarily any genetic data for those positions. So put in NA matrices. 
			if(positions == "all")
			{
				wgaimObject$geno[[chromosome]]$map <- wgaimObject$geno[[chromosome]]$dist <- currentChrProbabilitiesMap
				wgaimObject$geno[[chromosome]]$data <- matrix(nrow = nrow(transformed$finals), ncol = length(currentChrProbabilitiesMap), data = NA)
				rownames(wgaimObject$geno[[chromosome]]$data) <- rownames(transformed$finals)
				colnames(wgaimObject$geno[[chromosome]]$data) <- names(currentChrProbabilitiesMap)

				wgaimObject$geno[[chromosome]]$founders <- matrix(nrow = nrow(transformed$founders), ncol = length(currentChrProbabilitiesMap), data = NA)
				rownames(wgaimObject$geno[[chromosome]]$founders) <- rownames(transformed$founders)
				colnames(wgaimObject$geno[[chromosome]]$founders) <- names(currentChrProbabilitiesMap)

				markerIndices <- match(names(currentChrProbabilitiesMap), flattenedProbabilitiesMap)
				start <- min(markerIndices)
				end <- max(markerIndices)
				wgaimObject$geno[[chromosome]]$imputed.data <- transformed$probabilities[, (start*currentDataFounders - (currentDataFounders - 1)):(currentDataFounders*end)]
			}
			else
			{
				wgaimObject$geno[[chromosome]]$map <- wgaimObject$geno[[chromosome]]$dist <- currentChrMap
				wgaimObject$geno[[chromosome]]$data <- transformed$finals[,names(currentChrMap)]
				wgaimObject$geno[[chromosome]]$founders <- transformed$founders[,names(currentChrMap)]

				markerIndices <- match(names(currentChrMap), flattenedProbabilitiesMap)
				probabilityColumnIndices <- rep(markerIndices, times = currentDataFounders)*currentDataFounders + rep(((-(currentDataFounders - 1)):0), times = length(markerIndices))
				wgaimObject$geno[[chromosome]]$imputed.data <- transformed$probabilities[, probabilityColumnIndices]
			}
			class(wgaimObject$geno[[chromosome]]) <- "A"
		}
		class(wgaimObject) <- c("mpInterval", "cross", "interval")
		wgaimObject$gen.type <- "mpMarker"
		results[[length(results)+1]] <- wgaimObject
	}
	if(length(results) == 1) return(results[[1]])
	return(results)
}
