context("Test generateIntervalMidPoints")

test_that("Test generateIntervalMidPoints",
	{
		map <- sim.map(len = c(100, 100), n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		map[[1]] <- c("Extra1" = 0, map[[1]][1:5], "Extra2" = 40, map[[1]][6:11], "Extra3" = 100)
		map[[2]] <- c("Extra4" = 0, map[[2]][1:5], "Extra5" = 40, map[[2]][6:11], "Extra6" = 100)
		f2Pedigree <- f2Pedigree(10)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		newPositions <- generateIntervalMidPoints(cross)
		expect_identical(names(newPositions[[1]]), paste0("Chr1Interval", seq(5, 95, 10)))
		expect_identical(names(newPositions[[2]]), paste0("Chr2Interval", seq(5, 95, 10)))
	})
