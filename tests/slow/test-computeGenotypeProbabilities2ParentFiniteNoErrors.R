context("genotype probability computation, two parents, finite selfing, no errors")
test_that("Test zero generations of intercrossing, codominant markers, no errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_true(all(genotypesFromProbabilities == result@geneticData[[1]]@finals))
			expect_true(all(result@geneticData[[1]]@probabilities@data[1:30, 1:30] == 1 | result@geneticData[[1]]@probabilities@data[1:30, 1:30] == 0))
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 100)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test zero generations of intercrossing, dominant markers, no errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.9)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 500)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test non-zero generations of intercrossing, codominant markers, no errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is hetrozygous. 
			expect_true(all(genotypesFromProbabilities == result@geneticData[[1]]@finals))
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 100, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 100, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = 100, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = 100, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
test_that("Test non-zero generations of intercrossing, dominant markers, no errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.83)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
test_that("Test zero generations of intercrossing, codominant markers, no errors, with extra positions",
	{
		sampleSize <- 500
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"])))/sampleSize, 0.95)
			#In this particular scenario recombinations are marked as occuring between extra and D1M27, rather than the first pair. 
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"])))/sampleSize, 0.75)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = sampleSize)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test zero generations of intercrossing, dominant markers, no errors, with extra positions",
	{
		sampleSize <- 500
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"])))/sampleSize, 0.95)
			#In this particular scenario recombinations are marked as occuring between extra and D1M27, rather than the first pair. 
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"])))/sampleSize, 0.75)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = sampleSize)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test non-zero generations of intercrossing, codominant markers, no errors, with extra positions",
	{
		sampleSize <- 500
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"])))/sampleSize, 0.95)
			#In this particular scenario recombinations are marked as occuring between extra and D1M27, rather than the first pair. 
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"])))/sampleSize, 0.75)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
test_that("Test non-zero generations of intercrossing, dominant markers, no errors, with extra positions",
	{
		sampleSize <- 1000
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"])))/sampleSize, 0.75)
			#In this particular scenario recombinations are marked as occuring between extra and D1M27, rather than the first pair. 
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"])))/sampleSize, 0.745)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
