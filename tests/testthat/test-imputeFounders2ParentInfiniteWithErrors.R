context("Founder imputation, two parents, infinite selfing, with errors")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(map)
		{
			pedigree <- rilPedigree(500, selfingGenerations = 10)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0.05))

			#Imputed data may not be identical to orignal data, even though markers are fully informative
			tmp <- table(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)

			#Dominance doesn't really make a difference, because it's assumed inbred
			cross <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross, map = map)
			result <- imputeFounders(mapped, errorProb = 0.05)
			tmp <- table(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)
			
			errors <- result@geneticData[[1]]@imputed@errors
			expect_lt(sum(errors) / length(errors), 0.01)
		}
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
	})
test_that("Test non-zero generations of intercrossing",
	{
		testFunc <- function(map)
		{
			pedigree <- twoParentPedigree(initialPopulationSize = 1000, selfingGenerations = 10, intercrossingGenerations = 2, nSeeds = 1)
			pedigree@selfing <- "infinite"
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0.05))

			#Imputed data may not be identical to orignal data, even though markers are fully informative
			tmp <- table(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)
		
			#Dominance doesn't really make a difference, because it's assumed inbred
			cross <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross, map = map)
			result <- imputeFounders(mapped, errorProb = 0.05)
			tmp <- table(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)
			
			errors <- result@geneticData[[1]]@imputed@errors
			expect_lt(sum(errors) / length(errors), 0.01)
		}
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
	})	
test_that("Test removal of deliberate errors",
	{
		testFunc <- function(map, intercrossingGenerations)
		{
			pedigree <- twoParentPedigree(initialPopulationSize = 500, selfingGenerations = 10, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			pedigree@selfing <- "infinite"
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			#Add an error
			mapped@geneticData[[1]]@finals[,50] <- 1L
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0.05))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			naIndices <- result@geneticData[[1]]@finals == 3
			result@geneticData[[1]]@finals[naIndices] <- NA
			result@geneticData[[1]]@imputed@data[naIndices] <- NA
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)
			
			errors <- result@geneticData[[1]]@imputed@errors
			expect_lt(sum(errors[,-50]) / length(errors[,-50]), 0.01)
			expect_gt(sum(errors[,50]), 100)

			#Dominance doesn't really make a difference, because it's assumed inbred
			cross <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross, map = map)
			#Add an error
			mapped@geneticData[[1]]@finals[,50] <- 1L
			result <- imputeFounders(mapped, errorProb = 0.05)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)
			
			errors <- result@geneticData[[1]]@imputed@errors
			expect_lt(sum(errors[,-50]) / length(errors[,-50]), 0.01)
			expect_gt(sum(errors[,50]), 100)
		}
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map, 0)
		testFunc(map, 2)
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map, 0)
		testFunc(map, 2)
	})	
