context("Test addition of mpcross objects in estimating recombination fractions, for biparental designs")

test_that("Checking that an f2 design, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)
	})
test_that("Checking that a ril design, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		#With infinite generations of selfing assumed
		rilPedigree <- rilPedigree(500, selfingGenerations = 10)
		cross <- simulateMPCross(map=map, pedigree=rilPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		suppressWarnings(multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE))
		suppressWarnings(singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE))
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)

		#Now with only finite generations of selfing assumed.
		rilPedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=rilPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)

	})
test_that("Checking that a 4-way ril design with a single funnel, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		fourParentPedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		cross <- simulateMPCross(map=map, pedigree=fourParentPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		suppressWarnings(multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE))
		suppressWarnings(singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE))
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)
	})
test_that("Checking that a 4-way design without selfing and with a single funnel, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		fourParentPedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		fourParentPedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=fourParentPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)

	})
test_that("Checking that a 4-way ril design with random funnels and two generations of intercrossing, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		fourParentPedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 2)
		cross <- simulateMPCross(map=map, pedigree=fourParentPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		suppressWarnings(multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE))
		suppressWarnings(singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE))
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)
	})
test_that("Checking that a 4-way design without selfing, with random funnels and two generations of intercrossing, when split into 100 different datasets, gives the same RF estimates",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		fourParentPedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		fourParentPedigree@selfing <- "finite"
		cross <- simulateMPCross(map=map, pedigree=fourParentPedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:10)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		expect_identical(multipleDatasetsRf@rf@theta, singleDatasetRf@rf@theta)
		expect_equal(multipleDatasetsRf@rf@lod, singleDatasetRf@rf@lod)
		expect_equal(multipleDatasetsRf@rf@lkhd, singleDatasetRf@rf@lkhd)

	})
