#' @export
addExtraMarker <- function(mpcrossMapped, newMarker, useOnlyExtraImputationPoints = TRUE)
{
	if(!inherits(newMarker, "mpcross"))
	{
		stop("Input newMarker must be an object of class mpcross")
	}
	if(!inherits(mpcrossMapped, "mpcrossMapped"))
	{
		stop("Input mpcrossMapped must be an object of class mpcrossMapped")
	}
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Input mpcrossMapped must contain a single data set")
	}
	if(is.null(mpcrossMapped@geneticData[[1]]@imputed))
	{
		stop("Imputation data must be available")
	}
	marginalNewMarker <- table(finals(newMarker))
	nObservations <- sum(marginalNewMarker)
	if(useOnlyExtraImputationPoints)
	{
		imputationGridResults <- mpcrossMapped@geneticData[[1]]@imputed@data[,extraImputationPoints(mpcrossMapped)]
		chiSquared <- apply(imputationGridResults, 2, function(x)
		{
			x[x > 8] <- NA
			observed <- table(x, finals(newMarker))
			marginalImputed <- table(x[!is.na(finals(newMarker))])
			expected <- outer(marginalImputed, marginalNewMarker) / nObservations
			return(sum((observed - expected)^2 / expected))
		})
		bestLocation <- names(which.max(chiSquared))
		return(bestLocation)
	}
	else
	{
		stop("This code path is not implemented yet")
	}
}
