context("genotype probability computation, four parents, infinite selfing")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(4*x-3):(4*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is hetrozygous. 
			expect_true(all((genotypesFromProbabilities == result@geneticData[[1]]@finals) | is.na(result@geneticData[[1]]@finals)))
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 6, intercrossingGenerations = 1, nSeeds = 1)
		pedigrees[[2]] <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 6, intercrossingGenerations = 1, nSeeds = 1)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "infinite"
			testFunc(map, pedigree)
		}
	})
test_that("Test non-zero generations of intercrossing",
	{
		testFunc <- function(map, pedigree)
		{
			#First check that with fully informative markers we get back the original data. 
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(4*x-3):(4*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is hetrozygous. 
			expect_true(all((genotypesFromProbabilities == result@geneticData[[1]]@finals) | is.na(result@geneticData[[1]]@finals)))
		}
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigrees <- list()
		pedigrees[[1]] <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 1)
		pedigrees[[2]] <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 2)
		pedigrees[[3]] <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 1)
		pedigrees[[4]] <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 2)
		for(pedigree in pedigrees)
		{
			pedigree@selfing <- "infinite"
			testFunc(map, pedigree)
		}
	})
