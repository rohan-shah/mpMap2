#' @title Convert to mpMap format
#' @description Convert to the format used by the original mpMap package. 
#' @details Converts an \code{mpcross} object to the format used by the original mpMap, the predecessor of this package. It is unlikely that this function will ever need to be used. 
#' @param mpcross The object of class \code{mpcross} to convert. 
#' @return An object with structure compatible with the older mpMap package. 
#' @export
toMpMap <- function(mpcross)
{
	inheritsNewMpcrossArgument(mpcross)
	if(length(mpcross@geneticData[[1]]) > 1)
	{
		stop("Only objects with a single genetic dataset can be converted")
	}
	if(.hasSlot(mpcross, "rf") && !is.null(mpcross@rf))
	{
		warning("Recombination fraction data will not be retained in the converted object")
	}
	if(.hasSlot(mpcross, "lg") && !is.null(mpcross@lg))
	{
		warning("Linkage group data will not be retained in the converted object")
	}
	if(mpcross@geneticData[[1]]@pedigree@selfing != "infinite")
	{
		stop("Cannot convert an object with finite generations of selfing")
	}
	retVal <- list()
	retVal$founders <- founders(mpcross)
	retVal$finals <- finals(mpcross)
	pedigree <- mpcross@geneticData[[1]]@pedigree
	observed <- rep(0, length(pedigree@lineNames))
	observed[match(rownames(finals(mpcross)), pedigree@lineNames)] <- 1

	retVal$pedigree <- data.frame(id = 1:length(pedigree@lineNames), Male = pedigree@father, Female = pedigree@mother, Observed = observed)
	rownames(retVal$pedigree) <- pedigree@lineNames

	retVal$id <- match(rownames(finals(mpcross)), pedigree@lineNames)
	retVal$fid <- match(rownames(founders(mpcross)), pedigree@lineNames)

	if(.hasSlot(mpcross, "map") && !is.null(mpcross@map)) retVal$map <- mpcross@map
	class(retVal) <- "mpcross"

	if(nFounders(mpcross) == 8) attr(retVal, "type") <- "ri8self"
	else if(nFounders(mpcross) == 4) attr(retVal, "type") <- "ri4self"
	else stop("Number of founders must be 4 or 8 to convert to mpMap format")

	return(retVal)
}
