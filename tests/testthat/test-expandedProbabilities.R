context("Test two-locus probabilities")
	
test_that("Test that probabilities sum to 1",
{
	for(nFounders in c(2, 4, 8))
	{
		indices <- matrix(0, nrow = 0, ncol = 2)
		multiples <- c()
		for(i in 1:nFounders)
		{
			for(j in 1:nFounders)
			{
				for(k in 1:nFounders)
				{
					for(l in 1:k)
					{
						multiple <- 1
						if(k != l) multiple <- multiple * 2
						indices <- rbind(indices, c((i-1)*nFounders + j - 1 + 1, (k - 1)*nFounders + l - 1 + 1))
						multiples <- c(multiples, multiple)
					}
				}
			}
		}
		for(nFunnels in 1:2)
		{
			for(intercrossingGenerations in 0:1)
			{
				for(r in c(0, 0.25, 0.5))
				{
					for(selfingGenerations in c(0:1))
					{
						probabilities <- mpMap2:::expandedProbabilitiesFinite(nFounders = nFounders, nFunnels = nFunnels, intercrossingGenerations = intercrossingGenerations, selfingGenerations = selfingGenerations, r = r)
						expect_equal(sum(probabilities), 1)
						expect_true(isSymmetric(probabilities))
						#Check a specific symmetry - Essentially, if we're going to a heterozygous state, then we can count that state once, and double the probability. 
						expect_equal(sum(probabilities[indices]*multiples), 1)
					}
					probabilities <- mpMap2:::expandedProbabilitiesInfinite(nFounders = nFounders, nFunnels = nFunnels, intercrossingGenerations = intercrossingGenerations, r = r)
					expect_true(isSymmetric(probabilities))
					expect_equal(sum(probabilities), 1)					
				}
			}
		}
	}
})
