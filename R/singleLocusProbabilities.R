singleLocusProbabilities <- function(nFounders, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing)
{
	if(!(nFounders %in% c(2, 4, 8, 16))) stop("Input nFounders must be one of 2, 4, 8, or 16")
	if(infiniteSelfing) .Call("singleLocusProbabilitiesInfinite", nFounders, nFunnels, intercrossingGenerations, selfingGenerations, PACKAGE="mpMap2")
	else .Call("singleLocusProbabilitiesFinite", nFounders, nFunnels, intercrossingGenerations, selfingGenerations, PACKAGE="mpMap2")
}
