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
	if(is.logical(verbose))
	{
		if(is.na(verbose))
		{
			stop("Input verbose cannot be NA")
		}
		else if(verbose)
		{
			verbose <- list(verbose = TRUE, progressStyle = 3L)
		}
		else verbose <- list(verbose = FALSE, progressStyle = 3L)
	}
	else
	{
		if(!is.list(verbose) || !("progressStyle" %in% names(verbose)))
		{
			stop("Input verbose must be TRUE, FALSE, or a list with entries named progressStyle and verbose")
		}
		if(!("verbose" %in% names(verbose))) verbose$verbose <- TRUE
		if(length(verbose$progressStyle) != 1L || !(verbose$progressStyle %in% 0:3))
		{
			stop("Input verbose$progressStyle must have value 0, 1, 2 or 3")
		}
		if(!is.logical(verbose$verbose) || length(verbose$verbose) != 1L)
		{
			stop("Input verbose$verbose must have value 0, 1, 2 or 3")
		}
	}

	mpcrossLG@lg@imputedTheta <- list()
	for(counter in 1:length(mpcrossLG@lg@allGroups))
	{
		group <- mpcrossLG@lg@allGroups[counter]
		rawData <- .Call("imputeGroup", mpcrossLG, verbose, group)$theta
		mpcrossLG@lg@imputedTheta[[counter]] <- new("rawSymmetricMatrix", data = rawData, markers = names(which(mpcrossLG@lg@groups == group)), levels = mpcrossLG@rf@theta@levels)
	}
	names(mpcrossLG@lg@imputedTheta) <- as.character(mpcrossLG@lg@allGroups)
	return(mpcrossLG)
}
