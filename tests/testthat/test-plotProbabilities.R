context("Test plotProbabilities function")

test_that("Test that plotProbabilities gives informative errors",
{
	map <- qtl::sim.map(len = 5, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE, n.mar = 21)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "infinite"
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 2) + removeHets()
	crossMapped <- new("mpcrossMapped", cross, map = map)
	probabilities <- computeGenotypeProbabilities(crossMapped)

	expect_error(plotProbabilities(probabilities@geneticData[[1]]@probabilities), "Input object must have class mpcrossMapped or geneticData")
	expect_error(plotProbabilities(crossMapped@geneticData[[1]]), "No probabilities were found")

})
test_that("Test that plotProbabilities actually works",
{
	map <- qtl::sim.map(len = c(5, 5), anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE, n.mar = 21)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "infinite"
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 2) + removeHets()
	crossMapped <- new("mpcrossMapped", cross, map = map)
	probabilities <- computeGenotypeProbabilities(crossMapped, extraPositions = generateGridPositions(1.1))

	pdf(NULL)
		plotProbabilities(probabilities, chromosomes = "2")
	dev.off()
	pdf(NULL)
		expect_error(plotProbabilities(probabilities, positions = colnames(probabilities@geneticData[[1]]@probabilities@data)[8:9]), NA)
	dev.off()
})
