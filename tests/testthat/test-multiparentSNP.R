context("multiparentSNP")

test_that("Check that multiparentSNP works for a 4-way intercross",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + multiparentSNP(TRUE)
		cross3 <- cross + multiparentSNP(FALSE)
	})
test_that("Check that multiparentSNP doesn't work for an F2 or RIL",
	{
		pedigree <- f2Pedigree(1000)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(cross+ multiparentSNP(TRUE), throws_error())
		expect_that(cross+ multiparentSNP(FALSE), throws_error())
	})
