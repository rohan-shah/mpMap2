context("Test addition of mpcrossRf objects")

test_that("Cannot combine object with itself",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(10)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)
		expect_that(rf+rf, throws_error())
	})
test_that("Checking f2 pedigree split into 10 different datasets gives the same answer",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		separateRf <- estimateRF(subset(cross, lines = 1:50))
		combinedRf <- estimateRF(cross)
		for(i in 2:10)
		{
			suppressWarnings(capture.output(separateRf <- separateRf + estimateRF(subset(cross, lines = 1:50 + (i-1)*50))))
		}
		expect_identical(combinedRf@rf@theta, separateRf@rf@theta)
		expect_equal(combinedRf@rf@lod, separateRf@rf@lod, tolerance = 0.001)
		expect_equal(combinedRf@rf@lkhd, separateRf@rf@lkhd, tolerance = 0.001)
	})
test_that("Checking that f2 experiment split into two different subsets gives the same answer",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(1000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross@geneticData[[1]]@finals[1:500,1:3] <- cross@geneticData[[1]]@finals[501:1000,9:11] <- NA

		cross1 <- subset(cross, markers = 4:11)
		cross1 <- subset(cross1, lines = 1:500)

		cross2 <- subset(cross, markers = 1:8)
		cross2 <- subset(cross2, lines = 501:1000)

		rf1 <- estimateRF(cross1, keepLod = TRUE, keepLkhd = TRUE)
		rf2 <- estimateRF(cross2, keepLod = TRUE, keepLkhd = TRUE)

		#Two different recombination fraction calculations
		suppressWarnings(capture.output(combinedRf <- rf1 + rf2))
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		#After we give them the same marker orderings, they should have identical data, or close to it.
		combinedRf <- subset(combinedRf, markers = markers(rf))
		expect_identical(combinedRf@rf@theta@data, rf@rf@theta@data)
		expect_equal(combinedRf@rf@lod@x, rf@rf@lod@x)
		expect_equal(combinedRf@rf@lkhd@x, rf@rf@lkhd@x)
	})
