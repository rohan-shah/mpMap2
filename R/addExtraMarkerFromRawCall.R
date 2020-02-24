#' @include map-class.R
splitVector <- function(vector, splitValue)
{
	index <- match(splitValue, vector)
	if(index == length(vector)) return(list(before = vector[1:length(vector)], after = c()))
	return(list(before = vector[1:index], after = vector[(index+1):length(vector)]))
}
#' @title Add an extra marker from raw calling data
#' @description Add an extra marker to a map, based on raw calling data, using a QTL-mapping style approach. 
#' @param mpcrossMapped An object of class \code{mpcrossMapped} (dataset with a map), which must include imputed IBD genotypes and recombination fraction data. 
#' @param newMarker A matrix containing the raw data for the marker to add. 
#' @param useOnlyExtraImputationPoints Should we only attempt to add the new marker at points at which imputation data has been generated, which are \emph{not} markers? 
#' @return An object of class \code{addExtraMarkersStatistics} containing the test statistic values and the genetic map used to generate them. 
#' 
#' @details This function uses a QTL-mapping style approach to test for where an extra marker should be added to an existing map. The code uses the imputation data at a collection of points, and the \emph{raw calling data} for the extra marker. The raw calling data must be bivariate.
#'
#' Test statistics measuring the association of the new marker to a point are computed using a multivariate analysis of variance approach. If the imputed genotype at a point is independent of the data for the new marker, then the new marker probably should \emph{not} be mapped to that point. If the imputed genotype at a point and the data for the new marker are strongly \emph{dependent}, then the new marker \emph{should} probably be mapped to that point. Dependence and independence are measured using an F-test. 
#'
#' By default the set of points at which the new marker is considered for addition is the set of points at which imputation data is available, \emph{which are not markers}. The intention is that this set of points should be an equally spaced grid of points; this reduces the number of tests that are performed, as generally there are far fewer points in the grid, than there are markers. After the new marker is added, local reordering will need to be performed anyway, making any loss in accuracy by using the grid of points largely irrelevant. Setting \code{useOnlyExtraImputationPoints} to \code{FALSE} means that every marker position will also be used as a possible position for the new marker (this is not recommended). 
#' 
#' @export
addExtraMarkerFromRawCall <- function(mpcrossMapped, newMarker, useOnlyExtraImputationPoints = TRUE)
{
	if(!inherits(mpcrossMapped, "mpcrossMapped"))
	{
		stop("Input mpcrossMapped must be an object of class mpcrossMapped")
	}
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Input mpcrossMapped must contain a single data set")
	}
	if(any(!(rownames(newMarker) %in% lineNames(mpcrossMapped))))
	{
		stop("Invalid line names in new data")
	}
	if(is.null(mpcrossMapped@geneticData[[1]]@imputed))
	{
		stop("Imputation data must be available")
	}
	if(inherits(newMarker, "data.frame"))
	{
		stop("Input newMarker must be a matrix")
	}
	founders <- nFounders(mpcrossMapped)
	if(useOnlyExtraImputationPoints)
	{
		if(length(extraImputationPoints(mpcrossMapped)) == 0)
		{
			stop("There were no additional imputation points")
		}
		#Get out imputation results only for the grid points
		imputationGridResults <- imputationData(mpcrossMapped)[,extraImputationPoints(mpcrossMapped)]
		#Compute chi squared statistics
		testStatistics <- apply(imputationGridResults, 2, function(x)
		{
			model <- manova(newMarker ~ factor(x))
			return(summary(model)$stats[1, 3])
		})
		testMap <- imputationMap(mpcrossMapped)
		testMap <- lapply(testMap, function(x) x[names(x) %in% extraImputationPoints(mpcrossMapped)])
		class(testMap) <- "map"
		return(new("addExtraMarkersStatistics", data = testStatistics, map = testMap))
	}
	else
	{
		imputationResults <- imputationData(mpcrossMapped)[, flatImputationMapNames(mpcrossMapped)]
		#Compute chi squared statistics
		testStatistics <- apply(imputationResults, 2, function(x)
		{
			model <- manova(newMarker ~ factor(x))
			return(summary(model)$stats[1, 3])
		})
		testMap <- imputationMap(mpcrossMapped)
		return(new("addExtraMarkersStatistics", data = testStatistics, map = testMap))
		stop("This code path is not implemented yet")
	}
}
