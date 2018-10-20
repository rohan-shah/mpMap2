context("getChromosomes")

test_that("Checking that code to get chromosomes works",
	{
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		
		expect_identical(getChromosomes(mapped, c("D1M1", "D1M10", "D2M1", "D2M10")), c("D1M1" = "1", "D1M10" = "1", "D2M1" = "2", "D2M10" = "2"))
	})
