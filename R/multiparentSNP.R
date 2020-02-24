#' @include geneticData-class.R
#' @include mpcross-class.R
setClass("multiparentSNP", slots = list(keepHets = "logical"))
#' @title Convert all markers to SNP markers
#' @description Convert all markers in an object with fully informative markers to SNP markers
#' @details When initially generated, objects of class \code{mpcross} have markers that are fully informative - Every founder carries a different allele, and all marker heterozygotes are distinguishable. This function can be used to convert a simulated object to one with SNP markers. The resulting markers have two alleles, and the marker heterozygote might or might be observable.
#' @param keepHets Should heterozygotes for the SNP marker be kept?
#' @return An object of internal type \code{multiparentSNP}, which can be combined with an object of class \code{mpcross} using the addition operator. 
#' @export
multiparentSNP <- function(keepHets)
{
	return(new("multiparentSNP", keepHets = keepHets))
}
#' @rdname internalOperators
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
