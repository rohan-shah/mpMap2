splitVector <- function(vector, splitValue)
{
	index <- match(splitValue, vector)
	if(index == length(vector)) return(list(before = vector[1:length(vector)], after = c()))
	return(list(before = vector[1:index], after = vector[(index+1):length(vector)]))
}
#' @title Add an extra marker from raw calling data
#' @description Add an extra marker to a map, based on raw calling data, using a QTL-mapping style approach. 
#' @param mpcrossMapped An existing dataset with a map, which must include imputation data and recombination fraction data. 
#' @param newData A data.frame object containing the raw data for the marker to add
#' @param useOnlyExtraImputationPoints Use only the additional points at which imputation data has been generated? Currently this must be \code{TRUE}. 
#' @return A vector containing the test statistic values. 
#' 
#' @details This function uses a QTL-mapping style approach to test for where an extra marker should be added to an existing map. The code uses the imputation data at a collection of points, and the \emph{raw calling data} for the extra marker. Test statistics are computed using a multivariate analysis of variance approach. If the imputed genotype at a point is independent of the calling data, then that marker probably should not be mapped to that point. If the imputed genotype at a point and the calling data are strongly dependent, then the marker should probably be mapped to that point. Dependence and independence are measured using an F-test. 
#'
#' Currently the set of points which are tested is the set of points at which imputation data is available, \emph{which are not markers}. The intention is that this set of points should be an equally spaced grid of points; this has the affect of radically reducing the number of tests that are performed, as generally there are far fewer points in this grid, than there are markers. As the position chosen will need to be inspected and changed manually in any case, any loss in accuracy by using the grid of point is essentially irrelevant. In future it may be possible to use every marker position as the set of points at which tests are performed, by setting \code{useOnlyExtraImputationPoints} to \code{FALSE}. 
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
		stop("Input newMarker must be a data.frame")
	}
	founders <- nFounders(mpcrossMapped)
	if(useOnlyExtraImputationPoints)
	{
		if(length(extraImputationPoints(mpcrossMapped)) == 0)
		{
			stop("There were no additional imputation points")
		}
		#Get out the imputation map
		imputationMap <- imputationMap(mpcrossMapped)
		#Flatten the imputation map
		flattenedImputationMapNames <- flatImputationMapNames(mpcrossMapped)

		flattenedImputationMapPositions <- unlist(imputationMap)
		names(flattenedImputationMapPositions) <- flattenedImputationMapNames
		#Get out imputation results only for the grid points
		imputationGridResults <- imputationData(mpcrossMapped)[,extraImputationPoints(mpcrossMapped)]
		#Compute chi squared statistics
		testStatistics <- apply(imputationGridResults, 2, function(x)
		{
			model <- manova(newMarker ~ factor(x))
			return(summary(model)$stats[1, 3])
		})
		testMap <- imputationMap
		testMap <- lapply(testMap, function(x) x[names(x) %in% extraImputationPoints(mpcrossMapped)])
		class(testMap) <- "map"
		return(new("addExtraMarkersStatistics", data = testStatistics, map = testMap))
	}
	else
	{
		stop("This code path is not implemented yet")
	}
}
