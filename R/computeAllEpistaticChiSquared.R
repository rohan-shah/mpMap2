#' @title Compute chi-squared test statistics for independence
#' @description Compute chi-squared test statistics for independence
#' @details This function computes what are (approximately) chi-squared test statistics for independence of the genotypes at different points on the genome. This computation is done using the IBD probability data. Significant non-independence between IBD probabilities at distant points on the same chromosome, or points on different chromosomes, can indicate non-standard genetic inheritance or selective pressure. 
#' @param mpcrossMapped An object of class \code{mpcrossMapped} with IBD probability data. 
#' @param verbose If this is \code{TRUE} a progress bar is generated
#' @return A square matrix with rows and columns corresponding to genetic locations, and values corresponding to test statistics.
#' @export
computeAllEpistaticChiSquared <- function(mpcrossMapped, verbose = TRUE)
{
	if(length(mpcrossMapped@geneticData) > 1)
	{
		stop("Function computeAllEpistaticChiSquared can only be applied to a single design at a time")
	}
	if(is.null(mpcrossMapped@geneticData[[1]]@probabilities))
	{
		stop("Function computeAllEpistaticChiSquared requires probabilities to have been calculated")
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
	epistatic <- .Call("computeAllEpistaticChiSquared", mpcrossMapped@geneticData[[1]]@probabilities, nFounders(mpcrossMapped), selfing(mpcrossMapped@geneticData[[1]]@pedigree) == "infinite", verbose, PACKAGE="mpMap2")
	probabilitiesMap <- mpcrossMapped@geneticData[[1]]@probabilities@map
	rownames(epistatic) <- colnames(epistatic) <- unlist(lapply(probabilitiesMap, names))
	return(epistatic)
}
