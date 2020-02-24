condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call, ...)
  )
}
#' @title Impute missing recombination fraction estimates
#' @description Impute missing recombination fraction estimates
#' @details Recombination fractions between every pair of markers are estimated using numerical maximum likelihood. Unfortunately the likelihood is flat in some cases, so an estimate cannot be made. This later causes problems when trying to use estimated recombination fractions to order the markers, because a complete matrix of estimates is required. The solution is to impute the missing estimates using related estimates. For example, the recombination fraction between markers A and C may not be directly estimatable. However, there may be a marker B known to be tightly linked to A, which has a known recombination fraction with C. The estimated recombination fraction between B and C can be taken to be an estimate of the recombination fraction between A and C.
#' 
#' This function imputes values in the estimated recombination fraction matrix, to return a complete matrix. If there is a value that cannot be imputed, an error is triggered. Input \code{allErrors} controls whether the function will stop after encountering a single error, or continue and report all errors. If all errors are being reported, the optional function \code{extractErrorsFunction} is called with information about which missing estimates could not be imputed. 
#' @param mpcrossLG An object of class \code{mpcrossLG}, which contains estimated pairwise recombination fractions
#' @param verbose Should more verbose output be generated? 
#' @param allErrors If there is an error, should we immediately return, or should we continue, and report all errors?
#' @param extractErrorsFunction Error handling function. If there are errors and allErrors is \code{TRUE}, this function will be called with a matrix indicating which estimates could not be imputed. 
#' @return An object of class \code{mpcrossLG}, containing all the information in the input object, but also an imputed copy of the estimated recombination fraction data. 
#' @export
impute <- function(mpcrossLG, verbose = FALSE, allErrors = FALSE, extractErrorsFunction = function(e) e)
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
		withCallingHandlers({rawData <- .Call("imputeGroup", mpcrossLG, verbose, group, allErrors)$theta}, imputationErrors = extractErrorsFunction)
		mpcrossLG@lg@imputedTheta[[counter]] <- new("rawSymmetricMatrix", data = rawData, markers = names(which(mpcrossLG@lg@groups == group)), levels = mpcrossLG@rf@theta@levels)
	}
	names(mpcrossLG@lg@imputedTheta) <- as.character(mpcrossLG@lg@allGroups)
	return(mpcrossLG)
}
