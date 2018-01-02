context("Test single locus probabilities")
	
test_that("Test that probabilities sum to 1",
{
	for(nFounders in c(2, 4, 8))
	{
		indices <- matrix(0, nrow = 0, ncol = 2)
		multiples <- c()
		for(i in 1:nFounders)
		{
			for(j in 1:i)
			{
				multiple <- 1
				if(i != j) multiple <- multiple * 2
				indices <- rbind(indices, c(i, j))
				multiples <- c(multiples, multiple)
			}
		}
		for(selfingGenerations in c(0:1, Inf))
		{
			for(nFunnels in 1:2)
			{
				for(intercrossingGenerations in 0:1)
				{
					probabilities <- mpMap2:::singleLocusProbabilities(nFounders = nFounders, nFunnels = nFunnels, intercrossingGenerations = intercrossingGenerations, selfingGenerations = selfingGenerations, infiniteSelfing = selfingGenerations == Inf)
					#Check a specific symmetry - Essentially each het has two states, which contain half the probability. 
					expect_equal(sum(probabilities[indices] * multiples), 1)
					expect_equal(sum(probabilities), 1)
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
