context("Founder imputation, two parents, infinite selfing, no errors")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(map)
		{
			pedigree <- rilPedigree(500, selfingGenerations = 10)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			naIndices <- result@geneticData[[1]]@finals == 3
			result@geneticData[[1]]@finals[naIndices] <- NA
			result@geneticData[[1]]@imputed@data[naIndices] <- NA
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			#Dominance doesn't really make a difference, because it's assumed inbred
			cross <- cross + biparentalDominant()
			mapped <- new("mpcrossMapped", cross, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
		map <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		testFunc(map)
	})
