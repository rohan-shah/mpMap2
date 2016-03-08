context("Founder imputation, four parents, finite selfing")
test_that("Test zero generations of intercrossing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			expect_identical(result@geneticData[[1]]@imputed, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped)
			tmp <- table(result@geneticData[[1]]@imputed, cross@geneticData[[1]]@finals)
			endHet <- nrow(tmp)
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.8)
			#If a het is called, it should be correct with at least 68% chance
			expect_true(sum(diag(tmp[5:endHet,5:endHet])) / sum(tmp[5:endHet,]) > 0.68)
			#If a homozygote is called, it should be correct with a higher probability
			expect_true(sum(diag(tmp[1:4,1:4])) / sum(tmp[1:4,1:4]) > 0.95)
			#If a het is called, it should actually be a het (any het) with at least 75% chance
			expect_true(sum(tmp[5:endHet,5:endHet]) / sum(tmp[5:endHet,]) > 0.75)


			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped)
			tmp <- table(result@geneticData[[1]]@imputed, cross@geneticData[[1]]@finals)
			endHet <- nrow(tmp)
			#Correct imputation rate is kinda low
			expect_true(sum(diag(tmp)) / sum(tmp) > 0.5)
			#If a homozygote is called, it should be correct with a higher probability
			expect_true(sum(diag(tmp[1:4,1:4])) / sum(tmp[1:4,1:4]) > 0.95)
			#If a het is called, it should actually be a het (any het) with at least 70% chance
			if(nrow(tmp) < 5) browser()
			expect_true(sum(tmp[5:endHet,5:endHet]) / sum(tmp[5:endHet,]) > 0.70)
		}
		map1 <- sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "auto"
		pedigree2 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "auto"
		pedigree3 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
		pedigree3@selfing <- "auto"
		pedigree4 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
		pedigree4@selfing <- "auto"

		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}

	})
