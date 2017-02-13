context("changeMarkerPosition")

test_that("Checking that code to change marker positions works",
	{
		map <- qtl::sim.map(len = c(100, 100), n.mar = 11, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)

		altered <- changeMarkerPosition(mapped, marker = "D1M2", newChromosome = "2", newPosition = 0.5)
		expect_true(validObject(altered, complete=TRUE))
		expect_identical(markers(altered), c("D1M1", paste0("D1M", 3:11), "D2M1", "D1M2", paste0("D2M", 2:11)))

		altered <- changeMarkerPosition(mapped, marker = "D2M2", newChromosome = "1", newPosition = 0.5)
		expect_true(validObject(altered, complete=TRUE))
		expect_identical(markers(altered), c("D1M1", "D2M2", paste0("D1M", 2:11), "D2M1", paste0("D2M", 3:11)))
	})
