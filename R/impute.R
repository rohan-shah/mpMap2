#' @export
impute <- function(mpcrossLG, verbose = FALSE)
{
	isNewMpcrossLG(mpcrossLG)
	if(!is.null(mpcrossLG@lg@imputedTheta))
	{
		warning("Existing imputation data wil be overwritten")
	}
	mpcrossLG@lg@imputedTheta <- list()
	for(counter in 1:length(mpcrossLG@lg@allGroups))
	{
		group <- mpcrossLG@lg@allGroups[counter]
		mpcrossLG@lg@imputedTheta[[group]] <- .Call("imputeGroup", mpcrossLG, verbose, group)$theta
	}
	return(mpcrossLG)
}
