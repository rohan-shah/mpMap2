context("genotype probability computation, two parents, infinite selfing")
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
					apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is hetrozygous. 
			expect_true(all((genotypesFromProbabilities == result@geneticData[[1]]@finals) | is.na(result@geneticData[[1]]@finals)))
		}
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- rilPedigree(populationSize = 1000, selfingGenerations = 6)
		pedigree@selfing <- "infinite"
		testFunc(map, pedigree)
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
					apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			#The most probable founders should agree with the actual data, except for the case where the line really is hetrozygous. 
			expect_true(all((genotypesFromProbabilities == result@geneticData[[1]]@finals) | is.na(result@geneticData[[1]]@finals)))
		}
		map <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "infinite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "infinite"
		pedigrees <- list(pedigree1, pedigree2)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
