context("Test addition of mpcrossRF objects")

test_that("Cannot combine object with itself",
	{
		set.seed(1)
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(5000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		expect_that(cross+cross, throws_error())
	})
test_that("Checking f2 pedigree split into 100 different datasets",
	{
		set.seed(1)
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(5000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:100)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod, tolerance = 0.001)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd, tolerance = 0.001)
	})
