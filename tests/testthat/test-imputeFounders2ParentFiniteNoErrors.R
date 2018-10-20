context("Founder imputation, two parents, finite selfing, no errors")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(map, pedigree)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			#Dominance doesn't really make a difference, because it's assumed inbred
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.9)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)
		pedigree1 <- f2Pedigree(500)
		pedigree1@selfing <- "finite"
		pedigree2 <- rilPedigree(populationSize = 500, selfingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigrees <- list(pedigree1, pedigree2)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(map, pedigree)
			}
		}
	})
test_that("Test non-zero generations of intercrossing",
	{
		testFunc <- function(map, pedigree)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			#Dominance doesn't really make a difference, because it's assumed inbred
			cross2 <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.9)
			
			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)
		pedigree1 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "finite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 2, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "finite"
		pedigree3 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "finite"

		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(map, pedigree)
			}
		}
	})
