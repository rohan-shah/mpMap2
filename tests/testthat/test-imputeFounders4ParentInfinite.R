context("Founder imputation, four parents, infinite selfing")
test_that("Test zero generations of intercrossing, single chromosome",
	{
		testFunc <- function(pedigree)
		{
			map <- sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			naIndices <- result@geneticData[[1]]@finals > 4
			result@geneticData[[1]]@finals[naIndices] <- NA
			result@geneticData[[1]]@imputed[naIndices] <- NA
			expect_identical(result@geneticData[[1]]@imputed, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped)
			tmp <- table(result@geneticData[[1]]@imputed, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.96)
		}
		pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		testFunc(pedigree)
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		testFunc(pedigree)
	})
