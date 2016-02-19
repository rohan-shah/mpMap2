context("simulateMPCross")
test_that("Can simulate a cross object with a single markers",
	{
		pedigree <- f2Pedigree(100)
		map <- sim.map(len = 100, n.mar = 1, anchor.tel=FALSE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
	})
