context("Test constructLookupTable code")
test_that("Test two markers, two different marker patterns, no hets", 
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 6, intercrossingGenerations = 2)
	map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
	#Put in marker patterns for markers
	cross2 <- cross

	marker1Pattern <- c(1, 0, 0, 0, 0, 0, 0, 0)
	marker2Pattern <- c(0, 1, 0, 0, 0, 0, 0, 0)
	cross2@geneticData[[1]]@founders[,1] <- marker1Pattern
	cross2@geneticData[[1]]@founders[,2] <- marker2Pattern
	cross2@geneticData[[1]]@finals[,1] <- marker1Pattern[cross2@geneticData[[1]]@finals[,1]]
	cross2@geneticData[[1]]@finals[,2] <- marker2Pattern[cross2@geneticData[[1]]@finals[,2]]
	cross2@geneticData[[1]]@hetData[[1]] <- cross2@geneticData[[1]]@hetData[[2]] <- rbind(c(1, 1, 1), c(0, 0, 0))

	expect_identical(validObject(cross2, complete=TRUE), TRUE)
	recombValues <- c(0:20/200, 11:50/100)
	output <- capture.output(rf <- estimateRF(cross2, verbose = TRUE, recombValues = recombValues))
	#For example recombination value we have three tables, which are 2*2, of doubles
	expect_true(length(grep(paste0("Total lookup table size of ", length(recombValues)*3*2*2*8 , " bytes"), output)) == 1)
})
test_that("Test two markers, two different marker patterns, with hets", 
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 6, intercrossingGenerations = 2)
	map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
	#Put in marker patterns for markers
	cross2 <- cross

	marker1Pattern <- c(1, 0, 0, 0, 0, 0, 0, 0)
	marker2Pattern <- c(0, 1, 0, 0, 0, 0, 0, 0)
	cross2@geneticData[[1]]@founders[,1] <- marker1Pattern
	cross2@geneticData[[1]]@founders[,2] <- marker2Pattern
	cross2@geneticData[[1]]@finals[,1] <- marker1Pattern[cross2@geneticData[[1]]@finals[,1]]
	cross2@geneticData[[1]]@finals[,2] <- marker2Pattern[cross2@geneticData[[1]]@finals[,2]]
	cross2@geneticData[[1]]@hetData[[1]] <- cross2@geneticData[[1]]@hetData[[2]] <- rbind(c(1, 1, 1), c(0, 0, 0), c(1, 0, 2), c(0, 1, 2))

	expect_identical(validObject(cross2, complete=TRUE), TRUE)
	recombValues <- c(0:20/200, 11:50/100)
	output <- capture.output(rf <- estimateRF(cross2, verbose = TRUE, recombValues = recombValues))
	#For example recombination value we have three tables, which are 3*3, of doubles
	expect_true(length(grep(paste0("Total lookup table size of ", length(recombValues)*3*3*3*8 , " bytes"), output)) == 1)
})

test_that("Test two markers, single marker pattern, no hets", 
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 6, intercrossingGenerations = 2)
	map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
	#Put in marker patterns for markers
	cross2 <- cross

	marker1Pattern <- marker2Pattern <- c(1, 0, 0, 0, 0, 0, 0, 0)
	cross2@geneticData[[1]]@founders[,1] <- marker1Pattern
	cross2@geneticData[[1]]@founders[,2] <- marker2Pattern
	cross2@geneticData[[1]]@finals[,1] <- marker1Pattern[cross2@geneticData[[1]]@finals[,1]]
	cross2@geneticData[[1]]@finals[,2] <- marker2Pattern[cross2@geneticData[[1]]@finals[,2]]
	cross2@geneticData[[1]]@hetData[[1]] <- cross2@geneticData[[1]]@hetData[[2]] <- rbind(c(1, 1, 1), c(0, 0, 0))

	expect_identical(validObject(cross2, complete=TRUE), TRUE)
	recombValues <- c(0:20/200, 11:50/100)
	output <- capture.output(rf <- estimateRF(cross2, verbose = TRUE, recombValues = recombValues))
	#For example recombination value we have one table, 2*2, of doubles
	expect_true(length(grep(paste0("Total lookup table size of ", length(recombValues)*1*2*2*8 , " bytes"), output)) == 1)
})

test_that("Test two markers, single marker pattern, with hets", 
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 6, intercrossingGenerations = 2)
	map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
	#Put in marker patterns for markers
	cross2 <- cross

	marker1Pattern <- marker2Pattern <- c(1, 0, 0, 0, 0, 0, 0, 0)
	cross2@geneticData[[1]]@founders[,1] <- marker1Pattern
	cross2@geneticData[[1]]@founders[,2] <- marker2Pattern
	cross2@geneticData[[1]]@finals[,1] <- marker1Pattern[cross2@geneticData[[1]]@finals[,1]]
	cross2@geneticData[[1]]@finals[,2] <- marker2Pattern[cross2@geneticData[[1]]@finals[,2]]
	cross2@geneticData[[1]]@hetData[[1]] <- cross2@geneticData[[1]]@hetData[[2]] <- rbind(c(1, 1, 1), c(0, 0, 0), c(1, 0, 2), c(0, 1, 2))

	expect_identical(validObject(cross2, complete=TRUE), TRUE)
	recombValues <- c(0:20/200, 11:50/100)
	output <- capture.output(rf <- estimateRF(cross2, verbose = TRUE, recombValues = recombValues))
	#For example recombination value we have one table, 3*3, of doubles
	expect_true(length(grep(paste0("Total lookup table size of ", length(recombValues)*1*3*3*8 , " bytes"), output)) == 1)
})
