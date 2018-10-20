context("Founder imputation, eight parents, infinite selfing, no errors")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			naIndices <- result@geneticData[[1]]@finals > 8
			result@geneticData[[1]]@finals[naIndices] <- NA
			result@geneticData[[1]]@imputed@data[naIndices] <- NA
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.93)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 201, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 201, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- eightParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "infinite"
		pedigree2 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "infinite"
		pedigrees <- list(pedigree1, pedigree2)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}
	})
