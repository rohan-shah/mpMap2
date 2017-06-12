context("Test subsetting of mpcrossMapped objects")
test_that("Subsetting of mapped objects by lines discards rf data",
{
	map <- sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)
	expect_warning(subset(mapped, lines = rownames(finals(cross))[1:100]), "Subset function is discarding the recombination fraction data")

	subsettedMarkers <- subset(mapped, markers = names(mapped@map[[1]]))
	expect_is(subsettedMarkers, "mpcrossLG")
	expect_true(!is.null(subsettedMarkers@rf))
})
test_that("Subsetting of mapped objects by chromosomes keeps map and rf data",
{
	map <- sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)

	subsettedChromosomes <- subset(mapped, chromosomes = "1")
	expect_is(subsettedChromosomes, "mpcrossMapped")
	expect_true(!is.null(subsettedChromosomes@rf))
})
test_that("Subsetting of mapped objects by line names works for probabilities, with infinite selfing",
{
	map <- sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 2, selfingGenerations = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "infinite"
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)
	probabilities <- computeGenotypeProbabilities(mapped)

	subsettedFirst <- subset(probabilities, lines = finalNames(mapped)[1])
	subsettedFirst <- subset(probabilities, lines = 1)

	subsettedSecond <- subset(probabilities, lines = finalNames(mapped)[2])
	subsettedSecond <- subset(probabilities, lines = 2)
	
	expect_identical(subsettedFirst@geneticData[[1]]@probabilities@data, mapped@geneticData[[1]]@probabilities@data[1:4,])
	expect_identical(subsettedSecond@geneticData[[1]]@probabilities@data, mapped@geneticData[[1]]@probabilities@data[5:8,])
})
test_that("Subsetting of mapped objects by line names works for probabilities, with finite selfing",
{
	map <- sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 2, selfingGenerations = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "finite"
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)
	probabilities <- computeGenotypeProbabilities(mapped)

	subsettedFirst <- subset(probabilities, lines = finalNames(mapped)[1])
	subsettedFirst <- subset(probabilities, lines = 1)

	subsettedSecond <- subset(probabilities, lines = finalNames(mapped)[2])
	subsettedSecond <- subset(probabilities, lines = 2)
	
	expect_identical(subsettedFirst@geneticData[[1]]@probabilities@data, mapped@geneticData[[1]]@probabilities@data[1:10,])
	expect_identical(subsettedSecond@geneticData[[1]]@probabilities@data, mapped@geneticData[[1]]@probabilities@data[11:20,])
})
