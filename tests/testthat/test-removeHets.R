context("removeHets")

test_that("Check that removeHets works in F2",
	{
		pedigree <- f2Pedigree(100)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		withoutHets <- cross + removeHets()
		expect_equal(sum(is.na(withoutHets@geneticData[[1]]@finals)), sum(cross@geneticData[[1]]@finals == 3))
	})
