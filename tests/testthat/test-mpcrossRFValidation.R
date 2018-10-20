context("mpcrossRF validation")

test_that("Slot rf must have the right number of markers",
	{
		pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
		pedigree@selfing <- "finite"
		map1 <- qtl::sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		cross1 <- simulateMPCross(map = map1, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
		rf1 <- estimateRF(cross1)

		map2 <- qtl::sim.map(len = rep(100, 1), n.mar = 100, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		cross2 <- simulateMPCross(map = map2, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
		rf2 <- estimateRF(cross2)


		tmp1 <- rf1@rf
		tmp2 <- rf2@rf
		rf1@rf <- tmp2
		rf2@rf <- tmp1
		expect_that(validObject(rf1, complete=TRUE), throws_error())
		expect_that(validObject(rf2, complete=TRUE), throws_error())
	})
