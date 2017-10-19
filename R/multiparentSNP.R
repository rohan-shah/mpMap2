#' @include geneticData-class.R
#' @include mpcross-class.R
setClass("multiparentSNP", slots = list(keepHets = "logical"))
#' @export
multiparentSNP <- function(keepHets)
{
	return(new("multiparentSNP", keepHets = keepHets))
}
#' @rdname internalOperators
#' @description Internal operators, used to modify mpcross objects. 
setMethod(f = "+", signature = c("geneticData", "multiparentSNP"), definition = function(e1, e2)
{
	nFounders <- nFounders(e1)
	if(nFounders == 2)
	{
		stop("multiparentSNP cannot be applied to a biparental design. Did you mean to use biparentDominant?")
	}
	keepHets <- e2@keepHets
	if(keepHets)
	{
		return(multiparentSNPKeepHets(e1))
	}
	else return(multiparentSNPRemoveHets(e1))
})
multiparentSNPRemoveHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(apply(e1@founders, 2, function(x) length(unique(x))) != nFounders))
	{
		stop("Can only apply multiparentSNP to a fully informative multiparent design (one that has a different genotypes for every homozygote, at every marker)")
	}
	return(.Call("multiparentSNPRemoveHets", e1, PACKAGE="mpMap2"))

}
multiparentSNPKeepHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(unlist(lapply(e1@hetData, function(x) length(unique(x[,3])))) != nFounders + nFounders*(nFounders-1)/2))
	{
		stop("Can only apply multiparentSNP to a fully informative multiparent design (one that has nFounders*nFounders possible genotypes at every marker)")
	}
	return(.Call("multiparentSNPKeepHets", e1, PACKAGE="mpMap2"))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "multiparentSNP"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Assigning SNP marker patterns will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	if(any(nFounders(e1) == 2))
	{
		stop("multiparentSNP cannot be applied to a biparental design. Did you mean to use biparentDominant?")
	}
	e1@geneticData <- new("geneticDataList", lapply(e1@geneticData, function(x) x+e2))
	return(e1)
})
