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
	return(setdiff(unlist(lapply(mpcrossMapped@geneticData[[1]]@imputed@map, names)), unlist(lapply(mpcrossMapped@map, names))))
}
