context("Test warnImproperFunnels argument")
test_that("Test warnImproperFunnels with four-parent single funnel design",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)

		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), not(gives_warning()))

		pedigree@father[5] <- 1L
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), prints_text("Did you intend to use all"))

		pedigree@warnImproperFunnels <- FALSE
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), not(gives_warning()))

		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), not(gives_warning()))

		pedigree@father[5] <- 1L
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), throws_error("Repeated founders are only allowed with zero generations of intercrossing"))

		pedigree@father[5] <- 1L
		pedigree@warnImproperFunnels <- FALSE
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), throws_error("Repeated founders are only allowed with zero generations of intercrossing"))
	})
