expandedProbabilitiesInfinite <- function(nFounders, r, nFunnels, intercrossingGenerations)
{
	if(!(nFounders %in% c(2, 4, 8, 16))) stop("Input nFounders must be one of 2, 4, 8, or 16")
	.Call("expandedProbabilitiesInfinite", nFounders, r, nFunnels, intercrossingGenerations, PACKAGE="mpMap2")
}
expandedProbabilitiesFinite <- function(nFounders, r, nFunnels, intercrossingGenerations, selfingGenerations)
{
	if(!(nFounders %in% c(2, 4, 8, 16))) stop("Input nFounders must be one of 2, 4, 8, or 16")
	.Call("expandedProbabilitiesFinite", nFounders, r, nFunnels, intercrossingGenerations, selfingGenerations, PACKAGE="mpMap2")
}
