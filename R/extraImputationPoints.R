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
