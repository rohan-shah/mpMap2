context("getPositions")

test_that("Checking that code to get marker positions works",
	{
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		
		expect_identical(getPositions(mapped, c("D1M1", "D1M10", "D2M1", "D2M10")), c("D1M1" = 0, "D1M10" = 9, "D2M1" = 0, "D2M10" = 9))
	})
