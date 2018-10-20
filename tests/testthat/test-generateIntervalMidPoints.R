context("Test generateIntervalMidPoints")

test_that("Test generateIntervalMidPoints",
	{
		map <- qtl::sim.map(len = c(100, 100), n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		map[[1]] <- c("Extra1" = 0, map[[1]][1:5], "Extra2" = 40, map[[1]][6:11], "Extra3" = 100)
		map[[2]] <- c("Extra4" = 0, map[[2]][1:5], "Extra5" = 40, map[[2]][6:11], "Extra6" = 100)
		f2Pedigree <- f2Pedigree(10)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- new("mpcrossMapped", cross, map = map)
		newPositions <- generateIntervalMidPoints(cross)
		expect_identical(names(newPositions[[1]]), paste0("Chr1Interval", 1:length(newPositions[[1]])))
		expect_equal(newPositions[[1]], seq(5, 95, 10), check.attributes = FALSE)
		expect_identical(names(newPositions[[2]]), paste0("Chr2Interval", 1:length(newPositions[[2]])))
		expect_equal(newPositions[[2]], seq(5, 95, 10), check.attributes = FALSE)
	})
