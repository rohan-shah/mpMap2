#' @export 
combineKeepRF <- function(object1, object2, verbose = TRUE, gbLimit = -1, callEstimateRF = TRUE, skipValidity = FALSE)
{
	if(length(intersect(markers(object1), markers(object2))) != 0)
	{
		stop("Input objects to combineKeepRF must contain disjoint sets of markers")
	}
	if(is.null(object1@rf) || is.null(object2@rf))
	{
		stop("Input objects to combineKeepRF must contain recombination fraction data")
	}
	if(any(object1@rf@theta@levels != object2@rf@theta@levels))
	{
		stop("Slot rf@theta@levels must be the same for both input objects")
	}
	if(xor(is.null(object1@rf@lod), is.null(object2@rf@lod)) || xor(is.null(object1@rf@lkhd), is.null(object2@rf@lkhd)))
	{
		stop("Slots lod and lkhd must be present for both objects, or missing for both objects")
	}
	if(!is.logical(callEstimateRF) && !is.numeric(callEstimateRF))
	{
		stop("Input callEstimateRF must be TRUE, FALSE, or a numeric vector")
	}
	#Not worried about lineWeights warning
	suppressWarnings(combined <- addMpMap2(object1, object2, skipValidity = skipValidity))
	if(length(combined@geneticData) > 1)
	{
		stop("Could not combined objects into a single experiment")
	}
	newRF <- .Call("combineRFDisjoint", object1@rf, object2@rf, PACKAGE="mpMap2")
	combined <- new("mpcrossRF", combined, rf = newRF, skipValidity = skipValidity)
	combined@rf@gbLimit <- gbLimit
	#Not worried about warning for changing RF of existing object
	if(is.logical(callEstimateRF) && callEstimateRF)
	{
		suppressWarnings(combined <- estimateRF(combined, recombValues = object1@rf@theta@levels, markerRows = 1:nMarkers(combined), markerColumns = (nMarkers(object1)+1):nMarkers(combined), verbose = verbose, gbLimit = gbLimit))
	}
	else if(is.numeric(callEstimateRF))
	{
		suppressWarnings(combined <- estimateRF(combined, recombValues = object1@rf@theta@levels, markerRows = callEstimateRF, markerColumns = (nMarkers(object1)+1):nMarkers(combined), verbose = verbose, gbLimit = gbLimit))
	}
	return(combined)
}
