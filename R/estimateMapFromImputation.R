#' @title Re-estimate large gaps in a genetic map from IBD genotype imputations
#' @description Re-estimate large gaps in a genetic map from IBD genotype imputations
#' @details For larger gaps in a genetic map, the pairwise recombination fractions are not (by themselves) useful. An alternative is to estimate the IBD genotypes, and use the imputed IBD genotypes to re-estimate larger gaps using numerical maximum likelihood. Although the IBD genotypes are based on an existing genetic map, they may not be strongly affected by a large gap that has been poorly estimated, as the imputed IBD genotypes represent a consensus across all nearby markers, and also allow for genotyping errors. As a result, the re-estimated map may be different from the original map, and potentially more accurate. 
#' @param mpcrossMapped An object of class \code{mpcrossMapped}
#' @param gapSize The size of the gap to reestimate. 
#' @param recombinationFractions The recombination fractions to use for numerical maximum likelihood estimation
#' @return An object of class \code{mpcrossMapped} with a re-estimated map. 
#' @export
estimateMapFromImputation <- function(mpcrossMapped, gapSize = 5, recombinationFractions = c(0:60/600, 11:49/100))
{
	if(!inherits(mpcrossMapped, "mpcrossMapped"))
	{
		stop("Input object must have class mpcrossMapped")
	}
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Can only call estimateGapsInMap with a single experiment")
	}
	if(mpcrossMapped@geneticData[[1]]@pedigree@selfing == "finite")
	{
		stop("This function is only supported with infinite generations of selfing")
	}
	nFounders <- nFounders(mpcrossMapped)
	#The names of the markers at the right-hand endpoints of the gaps. 
	gaps <- unlist(lapply(mpcrossMapped@map, function(x) names(which(diff(x) > gapSize))))
	if(is.null(mpcrossMapped@geneticData[[1]]@imputed))
	{
		stop("Imputation data is required for this function")
	}
	intercrossingSelfingData <- getIntercrossingAndSelfingGenerations(mpcrossMapped)
	intercrossingGenerations <- unique(intercrossingSelfingData[,"intercrossing"])
	funnels <- getAllFunnels(mpcrossMapped, standardised = FALSE)
	#Used to determine whether we have only one funnel, or randomly chosen funnels. 
	nFunnels <- nrow(unique(funnels, MARGIN = 1))
	models <- lapply(as.list(intercrossingGenerations), 
		function(nIntercrossing)
		{
			sapply(recombinationFractions, 
				function(r)
				{
					expandedProbabilitiesInfinite(nFounders = nFounders(mpcrossMapped), r = r, nFunnels = nFunnels, intercrossingGenerations = nIntercrossing)
				}, simplify=FALSE)
		}
	)
	subsets <- sapply(intercrossingGenerations, function(x) subset(mpcrossMapped, lines = which(intercrossingSelfingData[,"intercrossing"] == x)), simplify = FALSE)
	for(rightHandMarker in gaps)
	{
		chromosome <- names(mpcrossMapped@map)[unlist(lapply(mpcrossMapped@map, function(x) rightHandMarker %in% names(x)))]
		rightHandMarkerIndex <- match(rightHandMarker, names(mpcrossMapped@map[[chromosome]]))
		leftHandMarker <- names(mpcrossMapped@map[[chromosome]])[rightHandMarkerIndex-1]
		jointLikelihood <- numeric(length(recombinationFractions))
		for(subsetIndex in 1:length(subsets))
		{
			leftHandData <- factor(imputationData(subsets[[subsetIndex]])[,leftHandMarker], levels = 1:nFounders)
			rightHandData <- factor(imputationData(subsets[[subsetIndex]])[,rightHandMarker], levels = 1:nFounders)
			jointTable <- table(leftHandData, rightHandData)
			jointLikelihood <- jointLikelihood + sapply(1:length(recombinationFractions), function(rIndex) sum(jointTable * log(models[[subsetIndex]][[rIndex]])))
		}
		newGapSize <- rfToHaldane(recombinationFractions[which.max(jointLikelihood)])
		mpcrossMapped@map[[chromosome]][rightHandMarkerIndex:length(mpcrossMapped@map[[chromosome]])] <- mpcrossMapped@map[[chromosome]][rightHandMarkerIndex:length(mpcrossMapped@map[[chromosome]])] - mpcrossMapped@map[[chromosome]][rightHandMarkerIndex] + newGapSize + mpcrossMapped@map[[chromosome]][rightHandMarkerIndex-1]
	}
	return(mpcrossMapped)
}
