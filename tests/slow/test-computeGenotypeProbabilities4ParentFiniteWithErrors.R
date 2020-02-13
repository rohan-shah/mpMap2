context("genotype probability computation, four parents, finite selfing, with errors")
test_that("Test fully informative markers, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:3, nLines(mapped@geneticData[[1]]), replace=TRUE)
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_gt(sum(diag(table(genotypesFromProbabilities, result@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.98)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA not assumed biased towards hets, with errors, no extra positions, no selfing",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:1, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.94)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA not assumed biased towards hets, with errors, no extra positions, with selfing",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:1, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.71)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA assumed biased towards hets, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:1, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, heterozygoteMissingProb = 0.9, homozygoteMissingProb = 0.01, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.77)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, with hets called, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.94)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, with hets called, with errors, no extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.72)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test fully informative markers, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(1:4, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.1))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nLines(result), 0.79)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M27"]))) / nLines(result), 0.79)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA not assumed biased towards hets, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:1, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nLines(result), 0.68)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M27"]))) / nLines(result), 0.68)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA assumed biased towards hets, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:1, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, heterozygoteMissingProb = 0.9, homozygoteMissingProb = 0.01, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.1))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nLines(result), 0.66)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M27"]))) / nLines(result), 0.66)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, with hets called, with errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- sample(0:2, nLines(mapped@geneticData[[1]]), replace=TRUE)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(10*x-9):(10*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M26"]))) / nLines(result), 0.68)
			expect_gt(sum(diag(table(genotypesFromProbabilities[,"extra"], cross@geneticData[[1]]@finals[,"D1M27"]))) / nLines(result), 0.68)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- fourParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- fourParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
