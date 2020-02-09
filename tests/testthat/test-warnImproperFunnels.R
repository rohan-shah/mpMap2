context("Test warnImproperFunnels argument")
test_that("Test warnImproperFunnels with four-parent single funnel design",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)

		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_warning(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), NA)

		pedigree@father[5] <- 1L
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), prints_text("Did you intend to use all"))

		pedigree@warnImproperFunnels <- FALSE
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_warning(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), NA)

		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_warning(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), NA)

		pedigree@father[5] <- 1L
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_error(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), "Repeated founders are only allowed with zero generations of intercrossing", class = "std::runtime_error")

		pedigree@father[5] <- 1L
		pedigree@warnImproperFunnels <- FALSE
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_error(rf <- estimateRF(cross, keepLod = FALSE, keepLkhd = FALSE), "Repeated founders are only allowed with zero generations of intercrossing", class = "std::runtime_error")
	})
