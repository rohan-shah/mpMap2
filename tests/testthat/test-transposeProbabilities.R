context("Test transposeProbabilities function")

test_that("Test transposeProbabilities with infinite selfing",
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 20, selfingGenerations = 1, intercrossingGenerations = 1)
	pedigree@selfing <- "infinite"
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(51, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
	mapped <- new("mpcrossMapped", grouped, map = map)
	probabilities <- computeGenotypeProbabilities(mapped)

	transposed <- transposeProbabilities(probabilities@geneticData[[1]])
	expect_identical(rownames(transposed), lineNames(probabilities))
	positionNames <- unlist(lapply(probabilities@geneticData[[1]]@probabilities@map, names))
	expect_identical(colnames(transposed), paste0(rep(positionNames, each = 8), " - ", rep(rownames(founders(probabilities)), times = length(positionNames))))
})
test_that("Test transposeProbabilities with finite selfing",
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 20, selfingGenerations = 1, intercrossingGenerations = 1)
	pedigree@selfing <- "finite"
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(51, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
	mapped <- new("mpcrossMapped", grouped, map = map)
	probabilities <- computeGenotypeProbabilities(mapped)

	transposed <- transposeProbabilities(probabilities@geneticData[[1]])
	expect_identical(rownames(transposed), lineNames(probabilities))
	positionNames <- unlist(lapply(probabilities@geneticData[[1]]@probabilities@map, names))
	expect_identical(colnames(transposed), paste0(rep(positionNames, each = 36), " - ", rep(1:36, times = length(positionNames))))
})
