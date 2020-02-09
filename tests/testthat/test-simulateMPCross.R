context("simulateMPCross")
test_that("Can simulate a cross object with a single markers",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 1, anchor.tel=FALSE, include.x=FALSE, eq.spacing=TRUE)
		expect_error(cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane), NA)
	})
test_that("Map is reordered",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		copiedMap <- map
		copiedMap[[1]] <- copiedMap[[1]][sample(1:101)]
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_identical(markers(cross), names(map[[1]]))
	})
