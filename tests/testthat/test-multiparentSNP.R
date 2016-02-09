context("multiparentSNP")
test_that("Check that multiparentSNP respects NA values",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 0)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		indicesNA <- sort(sample(1:(nMarkers(cross) * nLines(cross)), nMarkers(cross)*nLines(cross)/2, replace=FALSE))
		cross@geneticData[[1]]@finals[indicesNA] <- NA

		cross2 <- cross + multiparentSNP(keepHets = TRUE)
		expect_identical(which(is.na(cross2@geneticData[[1]]@finals)), indicesNA)

		cross@geneticData[[1]]@finals[cross@geneticData[[1]]@finals > 5] <- NA
		indicesNA <- which(is.na(cross@geneticData[[1]]@finals))
		cross3 <- cross + multiparentSNP(keepHets = FALSE)
		expect_identical(which(is.na(cross3@geneticData[[1]]@finals)), indicesNA)
	})

test_that("Check that multiparentSNP works for a 4-way intercross",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + multiparentSNP(TRUE)
		cross3 <- cross + multiparentSNP(FALSE)
	})
test_that("Check that multiparentSNP doesn't work for an F2 or RIL",
	{
		pedigree <- f2Pedigree(1000)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(cross+ multiparentSNP(keepHets = TRUE), throws_error())
		expect_that(cross+ multiparentSNP(keepHets = FALSE), throws_error())
	})
test_that("Check that multiparentSNP C code works the same as R code for four-way design", 
	{
		for(intecrossingGenerations in c(0, 3))
		{
			for(selfingGenerations in c(0, 6))
			{
				pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = selfingGenerations, nSeeds = 1, intercrossingGenerations = intecrossingGenerations)
				map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
				for(keepHets in c(FALSE, TRUE))
				{
					for(seed in 1:3)
					{
						cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
						set.seed(seed)
						snpMarkersC <- cross + multiparentSNP(keepHets = keepHets)
						set.seed(seed)
						snpMarkersR <- cross + mpMap2:::multiparentSNPPrototype(keepHets = keepHets)
						expect_identical(snpMarkersC, snpMarkersR)
					}
				}
			}
		}
	})
test_that("Check that multiparentSNP C code works the same as R code for eight-way design", 
	{
		for(intecrossingGenerations in c(0, 3))
		{
			for(selfingGenerations in c(0, 6))
			{
				pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = selfingGenerations, nSeeds = 1, intercrossingGenerations = intecrossingGenerations)
				map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
				for(keepHets in c(FALSE, TRUE))
				{
					for(seed in 1:3)
					{
						cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
						set.seed(seed)
						snpMarkersC <- cross + multiparentSNP(keepHets = keepHets)
						set.seed(seed)
						snpMarkersR <- cross + mpMap2:::multiparentSNPPrototype(keepHets = keepHets)
						expect_identical(snpMarkersC, snpMarkersR)
					}
				}
			}
		}
	})
