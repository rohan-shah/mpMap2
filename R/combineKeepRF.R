#' @export 
combineKeepRF <- function(object1, object2, verbose = TRUE, gbLimit = -1)
{
	if(length(intersect(markers(object1), markers(object2))) != 0)
	{
		stop("Input objects to combineKeepRF must contain disjoint sets of markers")
	}
	if(is.null(object1@rf) || is.null(object2@rf))
	{
		stop("Input objects to combineKeepRF must contain recombination fraction data")
	}
	if(xor(is.null(object1@rf@lod), is.null(object2@rf@lod)) || xor(is.null(object1@rf@lkhd), is.null(object2@rf@lkhd)))
	{
		stop("Slots lod and lkhd must be present for both objects, or missing for both objects")
	}
	#Not worried about lineWeights warning
	suppressWarnings(combined <- object1 + object2)
	if(length(combined@geneticData) > 1)
	{
		stop("Could not combined objects into a single experiment")
	}
	combined@rf <- .Call("combineRFDisjoint", object1@rf, object2@rf, PACKAGE="mpMap2")
	#Not worried about warning for changing RF of existing object
	suppressWarnings(combined <- estimateRF(combined, markerRows = 1:nMarkers(object1), markerColumns = (nMarkers(object1)+1):nMarkers(combined), verbose = verbose, gbLimit = gbLimit))
	return(combined)
}
