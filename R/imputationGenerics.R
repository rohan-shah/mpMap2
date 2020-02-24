#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Get names of positions for IBD genotype imputation
#' @rdname flatImputationMapNames
#' @description Get the names of all positions at which IBD genotype imputation has already been performed
#' @details Get the names of all positions at which IBD genotype imputation has already been performed
#' @param object The object from which to get the names of positions
#' @param ... Extra parameters, currently only \code{"experiment"} is supported. 
#' @return The names of all positions at which IBD genotype imputation has already been performed.
#' @export
setGeneric(name = "flatImputationMapNames", def = function(object, ...){standardGeneric("flatImputationMapNames")})
#' @rdname flatImputationMapNames
setMethod(f = "flatImputationMapNames", signature = "imputed", definition = function(object, ...)
{
	return(unlist(lapply(object@map, names)))
})
#' @rdname flatImputationMapNames
setMethod(f = "flatImputationMapNames", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Cannot extract the flattened imputation map names if there is no imputation data")
	}
	return(flatImputationMapNames(object@imputed))
})
#' @rdname flatImputationMapNames
setMethod(f = "flatImputationMapNames", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("experiment" %in% names(args) && args$experiment != 1) stop("Input experiment specified and not equal to one, but there is only one experiment in this object")
		return(flatImputationMapNames(object@geneticData[[1]]))
	}
	else
	{
		if("experiment" %in% names(args))
		{
			return(flatImputationMapNames(object@geneticData[[args$experiment]]))
		}
		else
		{
			return(lapply(object@geneticData, flatImputationMapNames))
		}
	}
})
#' @title Get map used for IBD genotype imputation
#' @rdname imputationMap
#' @description Get map used for IBD genotype imputation
#' @details Get the map of positions used for IBD genotype imputation. This is necessary because the points at which IBD genotype imputation has been performed may include non-marker points. See \code{\link{imputeFounders}} for further details. 
#' @param object The object from which to extract the IBD genotype imputation positions. 
#' @param ... Extra parameters. Currently only \code{"experiment"} is supported, letting the user extract the imputation map for a specific experiment. 
#' @return The map of positions used for IBD genotype imputation.
#' @export
setGeneric(name = "imputationMap", def = function(object, ...){standardGeneric("imputationMap")})
#' @rdname imputationMap
setMethod(f = "imputationMap", signature = "imputed", definition = function(object, ...)
{
	return(object@map)
})
#' @rdname imputationMap
setMethod(f = "imputationMap", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Cannot extract the imputation map if there is no imputation data")
	}
	return(imputationMap(object@imputed))
})
#' @rdname imputationMap
setMethod(f = "imputationMap", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("experiment" %in% names(args) && args$experiment != 1) stop("Input experiment specified and not equal to one, but there is only one experiment in this object")
		return(imputationMap(object@geneticData[[1]]))
	}
	else
	{
		if("experiment" %in% names(args))
		{
			return(imputationMap(object@geneticData[[args$experiment]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationMap))
		}
	}
})
#' @title Get out the IBD genotype imputation data
#' @description Get out the IBD genotype imputation data
#' @rdname imputationData
#' @details Extract the IBD genotype imputation data. The data takes the form of a matrix of values, with rows corresponding to genetic lines and columns corresponding to genetic positions. The genetic positions may include non-marker positions, so use \code{\link{imputationMap}} to find out the chromosome and position for every marker. 
#'
#' Each value in the matrix represents the predicted genotype for that genetic line, at that position. In the case of completely inbred experiments, each value in the matrix represents the founders from which that allele is believed to be derived. In the case of experiments with residual heterozygosity, the possible genotypes include heterozygotes, and the interpretation of the values in the matrix is more complicated. Function \code{\link{imputationKey}} gives information about how the values in the matrix correspond to actual genotypes. 
#' @param object The object from which to extract the IBD genotype imputation data
#' @param ... Extra parameters. Currently only \code{"experiment"} is supported, letting the user extract the imputation data for a specific experiment. 
#' @return The IBD genotype imputation data.
#' @export
setGeneric(name = "imputationData", def = function(object, ...){standardGeneric("imputationData")})
#' @rdname imputationData
setMethod(f = "imputationData", signature = "imputed", definition = function(object, ...)
{
	return(object@data)
})
#' @rdname imputationData
setMethod(f = "imputationData", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Imputation data is not present")
	}
	return(imputationData(object@imputed))
})
#' @rdname imputationData
setMethod(f = "imputationData", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("experiment" %in% names(args) && args$experiment != 1) stop("Input experiment specified and not equal to one, but there is only one experiment in this object")
		return(imputationData(object@geneticData[[1]]))
	}
	else
	{
		if("experiment" %in% names(args))
		{
			return(imputationData(object@geneticData[[args$experiment]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationData))
		}
	}
})
#' @title Get out key for IBD genotype imputations
#' @rdname imputationKey
#' @description Get out key for IBD genotype imputations
#' @details When IBD genotype imputation is performed using a population with finite generations of selfing, some of the imputed genotypes will be heterozygotes. However, the imputation code only returns a single value per line per genetic position. This key translates that value to a pair of founder alleles. 
#'
#' The key is a matrix with three columns. The first two columns represent founder alleles, and the third column gives the encoding for that particular pair of founder alleles.
#' @param object The object from which to get the imputation key. 
#' @param ... Extra parameters. Currently only \code{"experiment"} is supported, letting the user extract the imputation map for a specific experiment. 
#' @return Key giving the encoding of heterozygotes, in the imputed IBD genotype data. 
#' @examples
#' pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, 
#' 	selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 0)
#' selfing(pedigree) <- "finite"
#' #Generate map
#' map <- qtl::sim.map()
#' #Simulate data
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
#' crossSNP <- cross + multiparentSNP(keepHets = TRUE)
#' crossMapped <- mpcrossMapped(crossSNP, map = map)
#' imputed <- imputeFounders(crossMapped, errorProb = 0.01)
#' #An imputed IBD genotype of 1 indicates a homozygote for founder 1
#' #An imputed IBD genotype of 9 indicates a heterozygote for founders 1 and 2
#' #etc
#' head(imputationKey(imputed))
#' @export
setGeneric(name = "imputationKey", def = function(object, ...){standardGeneric("imputationKey")})
#' @rdname imputationKey
setMethod(f = "imputationKey", signature = "imputed", definition = function(object, ...)
{
	return(object@key)
})
#' @rdname imputationKey
setMethod(f = "imputationKey", signature = "geneticData", definition = function(object, ...)
{
	if(is.null(object@imputed))
	{
		stop("Imputation data is not present")
	}
	return(imputationKey(object@imputed))
})
#' @rdname imputationKey
setMethod(f = "imputationKey", signature = "mpcrossMapped", definition = function(object, ...)
{
	args <- list(...)
	if(length(object@geneticData) == 1)
	{
		if("experiment" %in% names(args) && args$experiment != 1) stop("Input experiment specified and not equal to one, but there is only one experiment in this object")
		return(imputationKey(object@geneticData[[1]]))
	}
	else
	{
		if("experiment" %in% names(args))
		{
			return(imputationKey(object@geneticData[[args$experiment]]))
		}
		else
		{
			return(lapply(object@geneticData, imputationKey))
		}
	}
})
