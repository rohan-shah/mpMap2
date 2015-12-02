context("Test addition of mpcrossRf objects")

test_that("Cannot combine object with itself",
	{
		set.seed(1)
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(5000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)
		expect_that(rf+rf, throws_error())
	})
test_that("Checking f2 pedigree split into 10 different datasets gives the same answer",
	{
		set.seed(1)
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		separateRf <- estimateRF(subset(cross, lines = 1:50))
		combinedRf <- estimateRF(cross)
		for(i in 2:10)
		{
			suppressWarnings(separateRf <- separateRf + estimateRF(subset(cross, lines = 1:50 + (i-1)*50)))
		}
		expect_identical(combinedRf@rf@theta, separateRf@rf@theta)
		expect_equal(combinedRf@rf@lod, separateRf@rf@lod, tolerance = 0.001)
		expect_equal(combinedRf@rf@lkhd, separateRf@rf@lkhd, tolerance = 0.001)
	})
