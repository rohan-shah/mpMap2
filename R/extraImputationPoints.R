#' @title Get out non-marker positions used for IBD genotype imputation
#' @description Get out non-marker positions used for IBD genotype imputation
#' @details Extract non-marker positions used for IBD genotype imputation
#' @param mpcrossMapped The object from which to get the non-marker positions
#' @return A vector of genetic position names. 
#' @export
extraImputationPoints <- function(mpcrossMapped)
{
	if(!inherits(mpcrossMapped, "mpcrossMapped"))
	{
		stop("Input object must have class mpcrossMapped")
	}
	if(length(mpcrossMapped) != 1)
	{
		stop("Input object must contain only a single experiment")
	}
	return(setdiff(flatImputationMapNames(mpcrossMapped), unlist(lapply(mpcrossMapped@map, names))))
}
