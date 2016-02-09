#' @export
impute <- function(mpcrossLG, verbose = FALSE)
{
	isNewMpcrossLGArgument(mpcrossLG)
	if(!is.null(mpcrossLG@lg@imputedTheta))
	{
		warning("Existing imputation data wil be overwritten")
	}
	if(is.null(mpcrossLG@rf))
	{
		stop("Input mpcrossLG object did not contain recombination fraction information")
	}
	mpcrossLG@lg@imputedTheta <- list()
	for(counter in 1:length(mpcrossLG@lg@allGroups))
	{
		group <- mpcrossLG@lg@allGroups[counter]
		rawData <- .Call("imputeGroup", mpcrossLG, verbose, group)$theta
		mpcrossLG@lg@imputedTheta[[counter]] <- new("rawSymmetricMatrix", data = rawData, markers = names(which(mpcrossLG@lg@groups == group)), levels = mpcrossLG@rf@theta@levels)
	}
	return(mpcrossLG)
}
