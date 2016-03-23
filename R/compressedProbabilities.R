compressedProbabilities <- function(nFounders, r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing)
{
	if(!(nFounders %in% c(2, 4, 8, 16))) stop("Input nFounders must be one of 2, 4, 8, or 16")
	.Call("compressedProbabilities", nFounders, r, nFunnels, intercrossingGenerations, selfingGenerations, infiniteSelfing, PACKAGE="mpMap2")
}
