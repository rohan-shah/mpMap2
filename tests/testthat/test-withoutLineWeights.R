context("Test that lineWeights argument works")
test_that("Checking that value of lineWeights option doesn't change results for f2",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(50)
		nReps <- 4
		for(i in 1:nReps)
		{
			cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
			allOneLineWeights <- estimateRF(subset(cross, lines = 1:49), keepLod = TRUE, keepLkhd = TRUE)
			differentLineWeights <- estimateRF(cross, lineWeights = c(rep(1, 49), 0), keepLod = TRUE, keepLkhd = TRUE)
			expect_identical(allOneLineWeights@rf@theta, differentLineWeights@rf@theta)
			expect_equal(allOneLineWeights@rf@lod, differentLineWeights@rf@lod)
			expect_equal(allOneLineWeights@rf@lkhd, differentLineWeights@rf@lkhd)
		}
	})
