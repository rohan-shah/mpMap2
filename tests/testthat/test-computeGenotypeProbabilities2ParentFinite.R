context("genotype probability computation, two parents, finite selfing")
test_that("Test zero generations of intercrossing, codominant markers",
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
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 1000)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test zero generations of intercrossing, dominant markers",
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
			expect_true(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities) > 0.85)
		}
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 1000)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test non-zero generations of intercrossing, codominant markers",
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
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
test_that("Test non-zero generations of intercrossing, dominant markers",
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
			expect_true(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities) > 0.85)
		}
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
