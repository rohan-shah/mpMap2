context("genotype probability computation, two parents, infinite selfing, no errors")
test_that("Test zero generations of intercrossing, no errors, no extra positions",
	{
		map <- sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
		}
		map <- sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
		map <- sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
		#The extra position should have essenitally thet same probabilities as the flanking markers
		expect_true(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman") > 0.92)
		expect_true(cor(genotypesFromProbabilities[,"extra"], genotypeFromProbabilities[,"D1M27"], method = "spearman") > 0.92)
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
			#The extra position should have essenitally thet same probabilities as the flanking markers
			expect_true(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman") > 0.92)
			expect_true(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M27"], method = "spearman") > 0.92)
		}
		map <- sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
