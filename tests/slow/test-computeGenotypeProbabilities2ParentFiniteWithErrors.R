context("genotype probability computation, two parents, finite selfing, with errors")
test_that("Test zero generations of intercrossing, codominant markers, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:3, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.5))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.9)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 500)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test zero generations of intercrossing, dominant markers, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.1))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			expect_gt(sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.73)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 1000)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test non-zero generations of intercrossing, codominant markers, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:3, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.9)
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
test_that("Test non-zero generations of intercrossing, dominant markers, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.67)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
test_that("Test zero generations of intercrossing, codominant markers, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:3, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.9)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 500)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test zero generations of intercrossing, dominant markers, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.1))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			#Our abilitiy to correctly impute hetrozygotes is very limited if there are no hets called!
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.5)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(populationSize = 1000)
		pedigree@selfing <- "finite"
		testFunc(map, pedigree)
	})
test_that("Test non-zero generations of intercrossing, codominant markers, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:3, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.9)
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
test_that("Test non-zero generations of intercrossing, dominant markers, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(3*x-2):(3*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)

			#Our abilitiy to correctly impute hetrozygotes is very limited if there are no hets called!		
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nrow(genotypesFromProbabilities), 0.5)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
