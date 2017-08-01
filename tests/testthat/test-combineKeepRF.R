context("combineKeepRF")

test_that("Test that combineKeepRF works",
{
	pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "finite"
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = "haldane", seed = 1) + multiparentSNP(keepHets = TRUE)
	rf <- estimateRF(cross)

	rf1 <- subset(rf, markers = 1:5)
	rf2 <- subset(rf, markers = 6:11)
	combined <- combineKeepRF(rf1, rf2, verbose = FALSE)
	expect_identical(rf, combined)

	rf1 <- subset(rf, markers = 1:10)
	rf2 <- subset(rf, markers = 11)
	combined <- combineKeepRF(rf1, rf2, verbose = FALSE)
	expect_identical(rf, combined)
})
