context("Test subsetting of mpcrossMapped objects")
test_that("Subsetting of mapped objects by markers discards map, unless keepMap is specified",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)
	
	suppressWarnings(subsetted <- subset(mapped, markers = 1:10))
	expect_equivalent(class(subsetted), "mpcrossLG")

	suppressWarnings(subsetted <- subset(mapped, markers = 1:10, keepMap = FALSE))
	expect_equivalent(class(subsetted), "mpcrossLG")

	suppressWarnings(subsetted <- subset(mapped, markers = 1:10, keepMap = TRUE))
	expect_equivalent(class(subsetted), "mpcrossMapped")
})
test_that("Subsetting of mapped objects by lines discards rf data",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
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
test_that("Subsetting of mapped objects by markers keeps rf data",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)

	subsettedChromosomes <- subset(mapped, markers = markers(mapped)[1:10])
	expect_is(subsettedChromosomes, "mpcrossLG")
	expect_true(!is.null(subsettedChromosomes@rf))
})
test_that("Subsetting of mapped objects by chromosomes keeps map and rf data",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
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
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 2, selfingGenerations = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "infinite"
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
	mapped <- new("mpcrossMapped", cross, map = map)
	probabilities <- computeGenotypeProbabilities(mapped)

	subsettedFirst <- subset(probabilities, lines = finalNames(mapped)[1])
	subsettedFirst <- subset(probabilities, lines = 1)

	subsettedSecond <- subset(probabilities, lines = finalNames(mapped)[2])
	subsettedSecond <- subset(probabilities, lines = 2)
	
	expect_identical(subsettedFirst@geneticData[[1]]@probabilities@data, probabilities@geneticData[[1]]@probabilities@data[1:4,])
	expect_identical(subsettedSecond@geneticData[[1]]@probabilities@data, probabilities@geneticData[[1]]@probabilities@data[5:8,])
})
test_that("Subsetting of mapped objects by line names works for probabilities, with finite selfing",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 2, selfingGenerations = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "finite"
	cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
	mapped <- new("mpcrossMapped", cross, map = map)
	probabilities <- computeGenotypeProbabilities(mapped)

	subsettedFirst <- subset(probabilities, lines = finalNames(mapped)[1])
	subsettedFirst <- subset(probabilities, lines = 1)

	subsettedSecond <- subset(probabilities, lines = finalNames(mapped)[2])
	subsettedSecond <- subset(probabilities, lines = 2)
	
	expect_identical(subsettedFirst@geneticData[[1]]@probabilities@data, probabilities@geneticData[[1]]@probabilities@data[1:10,])
	expect_identical(subsettedSecond@geneticData[[1]]@probabilities@data, probabilities@geneticData[[1]]@probabilities@data[11:20,])
})
