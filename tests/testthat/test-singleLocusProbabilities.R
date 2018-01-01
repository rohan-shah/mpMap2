context("Test addExtraMarkers")
	
test_that("Test that probabilities sum to 1",
{
	for(nFounders in c(2, 4, 8))
	{
		for(selfingGenerations in c(0:1, Inf))
		{
			for(nFunnels in 1:2)
			{
				for(intercrossingGenerations in 0:1)
				{
					expect_equal(sum(mpMap2:::singleLocusProbabilities(nFounders = nFounders, nFunnels = nFunnels, intercrossingGenerations = intercrossingGenerations, selfingGenerations = selfingGenerations, infiniteSelfing = selfingGenerations == Inf)), 1)
				}
			}
		}
	}
})
test_that("Test that F2 probabilities are correct", 
{
	probabilities <- mpMap2:::singleLocusProbabilities(nFounders = 2, nFunnels = 1, intercrossingGenerations = 0, selfingGenerations = 1, infiniteSelfing = FALSE)
	expected <- rep(0.25, 4)
	dim(expected) <- c(2, 2)
	expect_identical(probabilities, expected)
})
