context("genotype probability computation, two parents, infinite selfing, no errors")
test_that("Test zero generations of intercrossing, no errors, no extra positions",
	{
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- rilPedigree(populationSize = 1000, selfingGenerations = 6)
		pedigree@selfing <- "infinite"
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
		#Almost everything is a 1 or 0. The exception is hets, which end up coded as NA in the original dataset, and lead to probabilities that are neither 0 or 1. 
		booleans <- result@geneticData[[1]]@probabilities@data[1:100,1:20] == 1 | result@geneticData[[1]]@probabilities@data[1:100,1:20] == 0
		expect_gt(sum(booleans) / length(booleans), 0.92)
	})
test_that("Test non-zero generations of intercrossing, no errors, no extra positions",
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
			#Almost everything is a 1 or 0. The exception is hets, which end up coded as NA in the original dataset, and lead to probabilities that are neither 0 or 1. 
			booleans <- result@geneticData[[1]]@probabilities@data[1:100,1:20] == 1 | result@geneticData[[1]]@probabilities@data[1:100,1:20] == 0
			expect_gt(sum(booleans) / length(booleans), 0.92)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
test_that("Test zero generations of intercrossing, no errors, with extra positions",
	{
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- rilPedigree(populationSize = 1000, selfingGenerations = 6)
		pedigree@selfing <- "infinite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
		mapped <- new("mpcrossMapped", cross, map = map)
		suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
		genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
			{
				apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
			})
		genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))
		expect_true(all(result@geneticData[[1]]@probabilities@data[,"extra"] != 0) && all(result@geneticData[[1]]@probabilities@data[,"extra"] != 1))
		#Almost everything is a 1 or 0. The exception is hets, which end up coded as NA in the original dataset, and lead to probabilities that are neither 0 or 1. 
		booleans <- result@geneticData[[1]]@probabilities@data[1:100,1:20] == 1 | result@geneticData[[1]]@probabilities@data[1:100,1:20] == 0
		expect_gt(sum(booleans) / length(booleans), 0.92)
		#The extra position should have essenitally thet same probabilities as the flanking markers
		expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman"), 0.92)
		expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"], method = "spearman"), 0.92)
	})
test_that("Test non-zero generations of intercrossing, no errors, with extra positions",
	{
		testFunc <- function(map, pedigree)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5))))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))
			expect_true(all(result@geneticData[[1]]@probabilities@data[,"extra"] != 0) && all(result@geneticData[[1]]@probabilities@data[,"extra"] != 1))
			#Almost everything is a 1 or 0. The exception is hets, which end up coded as NA in the original dataset, and lead to probabilities that are neither 0 or 1. 
			booleans <- result@geneticData[[1]]@probabilities@data[1:100,1:20] == 1 | result@geneticData[[1]]@probabilities@data[1:100,1:20] == 0
			expect_gt(sum(booleans) / length(booleans), 0.85)
			#The extra position should have essenitally thet same probabilities as the flanking markers
			expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman"), 0.91)
			expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"], method = "spearman"), 0.91)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
