context("genotype probability computation, eight parents, finite selfing")
test_that("Test fully informative markers",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(36*x-35):(36*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is heterozygous. 
			expect_true(all(genotypesFromProbabilities == result@geneticData[[1]]@finals))
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA not assumed biased towards hets",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(36*x-35):(36*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_true(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities) > 0.68)
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[3]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, hets not called, and NA assumed biased towards hets",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets = FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, heterozygoteMissingProb = 0.9, homozygoteMissingProb = 0))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(36*x-35):(36*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.85)
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
test_that("Test SNP markers, with hets called",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(36*x-35):(36*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, most of the time. It won't be exact. 
			expect_gt(sum(diag(table(genotypesFromProbabilities, cross@geneticData[[1]]@finals))) / length(genotypesFromProbabilities), 0.85)
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[2]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[3]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[4]] <- eightParentPedigreeRandomFunnels(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[5]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[6]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[7]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 0, nSeeds = 1)
		pedigrees[[8]] <- eightParentPedigreeSingleFunnel(initialPopulationSize=500, selfingGenerations = 2, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "finite"
			testFunc(map, pedigree)
		}
	})
