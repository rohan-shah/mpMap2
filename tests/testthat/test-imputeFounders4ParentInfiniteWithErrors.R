context("Founder imputation, four parents, infinite selfing, with errors")
test_that("Test deliberate errors",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			errorIndices <- sample(1:length(mapped@geneticData[[1]]@finals), 0.01*length(mapped@geneticData[[1]]@finals))
			mapped@geneticData[[1]]@finals[errorIndices] <- sample(1:4, 0.01*length(mapped@geneticData[[1]]@finals), replace=TRUE) 
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0.01))
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.99)

			expect_gt(sum(result@geneticData[[1]]@imputed@errors), 0.5*0.01*length(mapped@geneticData[[1]]@finals))

			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			errorIndices <- sample(1:length(mapped@geneticData[[1]]@finals), 0.01 * length(mapped@geneticData[[1]]@finals))
			mapped@geneticData[[1]]@finals[errorIndices] <- sample(0:1, 0.01*length(mapped@geneticData[[1]]@finals), replace=TRUE)
			result <- imputeFounders(mapped, errorProb = 0.01)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.95)
			
			expect_gt(sum(result@geneticData[[1]]@imputed@errors), 0.25*0.01*length(mapped@geneticData[[1]]@finals))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "infinite"
		pedigree2 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "infinite"
		pedigree3 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 2)
		pedigree3@selfing <- "infinite"
		pedigree4 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 2)
		pedigree4@selfing <- "infinite"
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}
	})
