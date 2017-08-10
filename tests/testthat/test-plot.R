context("Test plot function")

test_that("Test that function plot gives informative errors",
{
	map <- qtl::sim.map(len = 5, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE, n.mar = 21)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "infinite"
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 2) + removeHets()
	crossMapped <- new("mpcrossMapped", cross, map = map)
	probabilities <- computeGenotypeProbabilities(crossMapped)
	imputed <- imputeFounders(crossMapped)

	expect_error(plot(crossMapped), "Function plot is not defined for an object of class \"mpcross\"")
	expect_error(plot(cross), "Function plot is not defined for an object of class \"mpcross\"")
	expect_error(plot(cross@geneticData[[1]]), "Function plot is not defined for an object of class \"geneticData\"")
	expect_error(plot(crossMapped@geneticData[[1]]), "Function plot is not defined for an object of class \"geneticData\"")
	expect_error(plot(probabilities@geneticData[[1]]@probabilities), "Function plot is not defined for an object of class \"probabilities\"")
	expect_error(plot(imputed@geneticData[[1]]@imputed), "Function plot is not defined for an object of class \"imputed\"")
})
